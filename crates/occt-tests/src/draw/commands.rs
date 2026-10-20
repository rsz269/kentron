//! DRAW command implementations, one method per command.

use kentron::geometry::curve::{BSplineCurve, BezierCurve, Circle, Curve, Ellipse, Hyperbola, Line, Parabola, TrimmedCurve};
use kentron::geometry::curve2::{BSplineCurve2, BezierCurve2, Circle2, Curve2, Ellipse2, Hyperbola2, Line2, Parabola2, TrimmedCurve2};
use kentron::geometry::surface::{BSplineSurface, BezierSurface, ExtrusionSurface, KnotVector, Plane, RevolutionSurface, Surface, TrimmedSurface};
use kentron::geometry::ParametricCurve;
use kentron::math::{Axis, Axis2, Direction, Direction2, Frame, Point, Point2, Transform, Transform2, Transformable, Vector, Vector2};
use kentron::primitives::{Cone, Cuboid, Cylinder, Sphere, Torus, Wedge};
use kentron::topology::{NodeKind, Orientation, Shape};
use kentron::sketch::{Outcome, Step};
use kentron::geometry::convert as gconvert;
use kentron::{blend, boolean, convert, make, modify, properties, sketch, sweep, tolerance, validity};

use super::{need, Draw, DrawError, DrawHistory, DrawObject, DrawResult};

/// A real as DRAW prints it (here the shortest form that reads back exactly).
fn format_real(x: f64) -> String {
    format!("{x}")
}

/// The arguments without the transformation commands' `-copy` and `-copymesh` options.
fn without_copy_flags<'a>(args: &[&'a str]) -> Vec<&'a str> {
    args.iter().copied().filter(|a| !matches!(*a, "-copy" | "-copymesh")).collect()
}

/// The transformation options of the arguments: `-copymesh` keeps the meshes.
fn transform_options(args: &[&str]) -> modify::TransformOptions {
    modify::TransformOptions { copy_meshes: args.contains(&"-copymesh") }
}

/// A real as DRAW prints tolerances in scientific notation: `1.00000e-07`.
fn format_c_exponent(x: f64) -> String {
    let text = format!("{x:.5e}");
    match text.split_once('e') {
        Some((mantissa, exponent)) => {
            let (sign, digits) = exponent.strip_prefix('-').map_or(("+", exponent), |d| ("-", d));
            format!("{mantissa}e{sign}{digits:0>2}")
        }
        None => text,
    }
}

// Shims appended by tools/draw2rs.py. Implement in place; the tool never
// modifies an existing method.
#[allow(non_snake_case)]
impl Draw {
    /// The degree, knots and multiplicities of a B-spline command's `name degree nknots k1 m1
    /// ...` arguments, and the numbers after them.
    fn bspline_knots(&self, args: &[&str], command: &str) -> DrawResult<(usize, Vec<f64>, Vec<usize>, Vec<f64>)> {
        need(args, 3, &format!("{command} name degree nknots knots... poles..."))?;
        let degree = self.number(args[1])? as usize;
        let nk = self.number(args[2])? as usize;
        let v = self.numbers(&args[3..])?;
        if v.len() < 2 * nk {
            return Err(DrawError::Failed(format!("{command}: missing knots")));
        }
        let knots = (0..nk).map(|i| v[2 * i]).collect();
        let mults = (0..nk).map(|i| v[2 * i + 1] as usize).collect();
        Ok((degree, knots, mults, v[2 * nk..].to_vec()))
    }
    /// DRAW command `2dapprox`: 1 uses in 1 converted scripts.
    ///
    /// Example: `2dapprox c 5 0 0 3 4 -1 4 -4 0 -4 -3`
    pub fn _2dapprox(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dapprox", args)
    }
    /// DRAW command `2dbeziercurve`: 13 uses in 10 converted scripts.
    ///
    /// Example: `2dbeziercurve c1 3 35 0 1 39 0 1 39 4 1`
    pub fn _2dbeziercurve(&mut self, args: &[&str]) -> DrawResult {
        // `2dbeziercurve name n x1 y1 [w1] ... xn yn [wn]`.
        need(args, 2, "2dbeziercurve name n poles...")?;
        let n = self.number(args[1])? as usize;
        let v = self.numbers(&args[2..])?;
        let stride = if v.len() == 3 * n { 3 } else if v.len() == 2 * n { 2 } else {
            return Err(DrawError::Failed("2dbeziercurve: wrong number of coordinates".into()));
        };
        let poles: Vec<Point2> = v.chunks(stride).map(|c| Point2::new(c[0], c[1])).collect();
        let curve = match stride {
            3 => BezierCurve2::rational(poles, v.chunks(3).map(|c| c[2]).collect()),
            _ => BezierCurve2::new(poles),
        };
        self.set(args[0], DrawObject::Curve2(Curve2::Bezier(curve)));
        Ok(String::new())
    }
    /// DRAW command `2dbsplinecurve`: 19 uses in 19 converted scripts.
    ///
    /// Example: `2dbsplinecurve cc 3 2 0 4 1 4 -1 -1 1 0 -1 1 0 0 1 0 0 1`
    pub fn _2dbsplinecurve(&mut self, args: &[&str]) -> DrawResult {
        // `2dbsplinecurve name degree nknots k1 m1 ... x1 y1 w1 ...`.
        let (degree, knots, mults, rest) = self.bspline_knots(args, "2dbsplinecurve")?;
        if rest.len() % 3 != 0 {
            return Err(DrawError::Failed("2dbsplinecurve: poles come as x y w".into()));
        }
        let poles: Vec<Point2> = rest.chunks(3).map(|c| Point2::new(c[0], c[1])).collect();
        let weights: Vec<f64> = rest.chunks(3).map(|c| c[2]).collect();
        let mut curve = BSplineCurve2::try_new(poles, knots, mults, degree)?;
        if weights.iter().any(|&w| w != 1.0) {
            curve = curve.with_weights(weights);
        }
        self.set(args[0], DrawObject::Curve2(Curve2::BSpline(curve)));
        Ok(String::new())
    }
    /// DRAW command `2dcvalue`: 13 uses in 11 converted scripts.
    ///
    /// Example: `2dcvalue c2d2 u1 x y1`
    pub fn _2dcvalue(&mut self, args: &[&str]) -> DrawResult {
        // `2dcvalue curve u x y [dx dy [d2x d2y]]`.
        need(args, 4, "2dcvalue curve u x y [dx dy [d2x d2y]]")?;
        let curve = match self.get(args[0]) {
            Some(DrawObject::Curve2(c)) => c.clone(),
            _ => return Err(DrawError::Failed(format!("{} is not a 2D curve", args[0]))),
        };
        let u = self.number(args[1])?;
        let (p, d1, d2) = curve.d2(u);
        let values = [p.x, p.y, d1.x, d1.y, d2.x, d2.y];
        for (name, value) in args[2..].iter().zip(values) {
            self.set(name, DrawObject::Number(value));
        }
        Ok(String::new())
    }
    /// DRAW command `2dextrema`: 1 uses in 1 converted scripts.
    ///
    /// Example: `2dextrema c2d e2d`
    pub fn _2dextrema(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dextrema", args)
    }
    /// DRAW command `2dinterpole`: 1 uses in 1 converted scripts.
    ///
    /// Example: `2dinterpole c 5 0 0 3 4 -1 4 -4 0 -4 -3`
    pub fn _2dinterpole(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dinterpole", args)
    }
    /// DRAW command `2dintersect`: 2 uses in 2 converted scripts.
    ///
    /// Example: `2dintersect l1 c1`
    pub fn _2dintersect(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dintersect", args)
    }
    /// DRAW command `2dlmirror`: 7 uses in 7 converted scripts.
    ///
    /// Example: `2dlmirror result 0 0 0 1`
    pub fn _2dlmirror(&mut self, args: &[&str]) -> DrawResult {
        // `2dlmirror name ... x y dx dy`.
        need(args, 5, "2dlmirror name ... x y dx dy")?;
        let n = args.len();
        let v = self.numbers(&args[n - 4..])?;
        let d = Direction2::try_new(v[2], v[3]).ok_or_else(|| DrawError::Failed("2dlmirror: null direction".into()))?;
        self.transform_vars2(&args[..n - 4], &Transform2::mirror_axis(Axis2::new(Point2::new(v[0], v[1]), d)))
    }
    /// DRAW command `2dpl`: 10 uses in 5 converted scripts.
    ///
    /// Example: `2dpl`
    pub fn _2dpl(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dpl", args)
    }
    /// DRAW command `2dpmirror`: 6 uses in 6 converted scripts.
    ///
    /// Example: `2dpmirror result 0 0`
    pub fn _2dpmirror(&mut self, args: &[&str]) -> DrawResult {
        // `2dpmirror name ... x y`.
        need(args, 3, "2dpmirror name ... x y")?;
        let n = args.len();
        let v = self.numbers(&args[n - 2..])?;
        self.transform_vars2(&args[..n - 2], &Transform2::mirror_point(Point2::new(v[0], v[1])))
    }
    /// DRAW command `2dpu`: 10 uses in 5 converted scripts.
    ///
    /// Example: `2dpu`
    pub fn _2dpu(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dpu", args)
    }
    /// DRAW command `2drotate`: 6 uses in 6 converted scripts.
    ///
    /// Example: `2drotate result 0 0 90`
    pub fn _2drotate(&mut self, args: &[&str]) -> DrawResult {
        // `2drotate name ... x y angle` (the angle in degrees).
        need(args, 4, "2drotate name ... x y angle")?;
        let n = args.len();
        let v = self.numbers(&args[n - 3..])?;
        self.transform_vars2(&args[..n - 3], &Transform2::rotation(Point2::new(v[0], v[1]), v[2].to_radians()))
    }
    /// DRAW command `2dtranslate`: 16 uses in 16 converted scripts.
    ///
    /// Example: `2dtranslate result 0 -4`
    pub fn _2dtranslate(&mut self, args: &[&str]) -> DrawResult {
        // `2dtranslate name ... dx dy`.
        need(args, 3, "2dtranslate name ... dx dy")?;
        let n = args.len();
        let v = self.numbers(&args[n - 2..])?;
        self.transform_vars2(&args[..n - 2], &Transform2::translation(Vector2::new(v[0], v[1])))
    }
    /// DRAW command `2dzoom`: 1 uses in 1 converted scripts.
    ///
    /// Example: `2dzoom 120`
    pub fn _2dzoom(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("2dzoom", args)
    }
    /// DRAW command `BUC60652`: 5 uses in 5 converted scripts.
    ///
    /// Example: `BUC60652 face`
    pub fn BUC60652(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("BUC60652", args)
    }
    /// DRAW command `BUC60792`: 1 uses in 1 converted scripts.
    ///
    /// Example: `BUC60792`
    pub fn BUC60792(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("BUC60792", args)
    }
    /// DRAW command `BUC60811`: 1 uses in 1 converted scripts.
    ///
    /// Example: `BUC60811`
    pub fn BUC60811(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("BUC60811", args)
    }
    /// DRAW command `BUC60843`: 2 uses in 2 converted scripts.
    ///
    /// Example: `BUC60843 result test_circle test_ellipse 0.1`
    pub fn BUC60843(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("BUC60843", args)
    }
    /// DRAW command `BUC60944`: 1 uses in 1 converted scripts.
    ///
    /// Example: `BUC60944 ../../..`
    pub fn BUC60944(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("BUC60944", args)
    }
    /// DRAW command `DT_ShapeDivide`: 4 uses in 4 converted scripts.
    ///
    /// Example: `DT_ShapeDivide ts ts0`
    pub fn DT_ShapeDivide(&mut self, args: &[&str]) -> DrawResult {
        // `DT_ShapeDivide result shape [tolerance]`: a copy divided where its curves and
        // surfaces are not C1 (the tolerance not taken into account).
        if !(2..=3).contains(&args.len()) {
            return self.not_implemented("DT_ShapeDivide", args);
        }
        let shape = self.shape_var(args[1])?;
        let copy = kentron::modify::copy(&mut self.model, shape, false).unwrap_or(shape);
        let divided = kentron::heal::split_at_breaks(&mut self.model, copy)?;
        self.set_shape_var(args[0], divided)?;
        Ok("Status: OK\n".into())
    }
    /// DRAW command `DT_ShapeConvertRev`: OCCT's turns the surfaces into surfaces of
    /// revolution and the curves and surfaces into Bézier pieces; kentron's converts them all
    /// to B-splines (its NURBS conversion), not cut into Bézier pieces.
    ///
    /// Example: `DT_ShapeConvertRev result a 1 1`
    pub fn DT_ShapeConvertRev(&mut self, args: &[&str]) -> DrawResult {
        // `DT_ShapeConvertRev result shape convert-2d convert-3d` (both 1 only).
        if args.len() != 4 || self.number(args[2])? != 1.0 || self.number(args[3])? != 1.0 {
            return self.not_implemented("DT_ShapeConvertRev", args);
        }
        let shape = self.shape_var(args[1])?;
        let converted = kentron::convert::to_nurbs(&mut self.model, shape)?;
        self.set_shape_var(args[0], converted)?;
        Ok("Status: DONE1\n".into())
    }
    /// DRAW command `expshape`: how many of a shape's surfaces, curves and pcurves are of
    /// each kind, in OCCT's report.
    ///
    /// Example: `expshape result 0 0`
    pub fn expshape(&mut self, args: &[&str]) -> DrawResult {
        // `expshape shape [degree [spans]]`: B-splines and Béziers of more than `degree`
        // (and B-splines of more than `spans` spans) counted, and those less than C1 (by
        // their representation).
        if args.is_empty() || args.len() > 3 {
            return Err(DrawError::Failed("usage: expshape shape [degree [spans]]".into()));
        }
        use kentron::geometry::surface::Surface;
        let shape = self.shape_var(args[0])?;
        let degree = match args.get(1) { Some(a) => self.number(a)? as usize, None => 0 };
        let spans = match args.get(2) { Some(a) => self.number(a)? as usize, None => 0 };
        let spans_of = |knots: &[f64]| knots.len().saturating_sub(1);
        let mut n = [0usize; 26];
        let faces = self.model.unique_subshapes(shape, NodeKind::Face);
        for &face in &faces {
            // Seams: the uses of an edge the face uses twice.
            let coedges = self.model.unique_subshapes(face, NodeKind::CoEdge);
            let edges: Vec<_> = coedges.iter().map(|c| self.model[kentron::topology::CoEdgeId::try_from(c.node).expect("a coedge")].edge).collect();
            n[0] += edges.iter().filter(|e| edges.iter().filter(|f| f == e).count() > 1).count();
            let Some(mut surface) = self.model.surface(face) else { continue };
            while let Surface::Trimmed(t) = surface {
                surface = t.basis();
            }
            match surface {
                Surface::BSpline(b) => {
                    n[1] += (b.u_degree().max(b.v_degree()) > degree) as usize;
                    n[2] += (spans_of(b.u_knots()).max(spans_of(b.v_knots())) > spans) as usize;
                    n[3] += (b.is_u_rational() || b.is_v_rational()) as usize;
                    let (b0, b1, b2, b3) = kentron::geometry::ParametricSurface::bounds(surface);
                    let (us, vs) = kentron::heal::surface_breaks(surface, b0, b1, b2, b3);
                    n[4] += (!us.is_empty() || !vs.is_empty()) as usize;
                }
                Surface::Bezier(b) => {
                    n[5] += (b.u_degree().max(b.v_degree()) > degree) as usize;
                    n[6] += (b.is_u_rational() || b.is_v_rational()) as usize;
                }
                Surface::Plane(_) => n[8] += 1,
                _ => n[9] += 1,
            }
            for c in &coedges {
                let Some((mut p, _, _)) = self.model.effective_pcurve(kentron::topology::CoEdgeId::try_from(c.node).expect("a coedge")) else { continue };
                while let Curve2::Trimmed(t) = p {
                    p = t.basis().clone();
                }
                match &p {
                    Curve2::BSpline(b) => {
                        n[18] += (b.degree() > degree) as usize;
                        n[19] += (spans_of(b.knots()) > spans) as usize;
                        n[20] += b.is_rational() as usize;
                        let k = b.knots();
                        n[21] += (!kentron::heal::curve2_breaks(&p, k[0], k[k.len() - 1]).is_empty()) as usize;
                    }
                    Curve2::Bezier(b) => {
                        n[22] += (b.degree() > degree) as usize;
                        n[23] += b.is_rational() as usize;
                    }
                    _ => n[25] += 1,
                }
            }
        }
        for edge in self.model.unique_subshapes(shape, NodeKind::Edge) {
            let Some((mut c, _, _)) = self.model.edge_curve(edge) else { continue };
            while let Curve::Trimmed(t) = c {
                c = t.basis().clone();
            }
            match &c {
                Curve::BSpline(b) => {
                    n[10] += (b.degree() > degree) as usize;
                    n[11] += (spans_of(b.knots()) > spans) as usize;
                    n[12] += b.is_rational() as usize;
                    let k = b.knots();
                    n[13] += (!kentron::heal::curve_breaks(&c, k[0], k[k.len() - 1]).is_empty()) as usize;
                }
                Curve::Bezier(b) => {
                    n[14] += (b.degree() > degree) as usize;
                    n[15] += b.is_rational() as usize;
                }
                _ => n[17] += 1,
            }
        }
        let labels = [
            "Number of seam edges - ",
            "Number of BSpline surfaces with degree more then 0 - ",
            "Number of BSpline surfaces with number of spans more then 0 - ",
            "Number of Rational BSpline surfaces ",
            "Number of BSpline surfaces with continuity less than specified - ",
            "Number of Bezier surfaces with degree more then 0 - ",
            "Number of Rational Bezier surfaces  - ",
            "Number of Bezier surfaces with continuity less than specified - ",
            "Number of Planes - ",
            "Number of other surfaces - ",
            "Number of BSpline curves with degree more then 0 - ",
            "Number of BSpline curves with number of spans more then - 0 - ",
            "Number of Rational BSpline curves ",
            "Number of  BSpline curves with less continuity - ",
            "Number of Bezier curves with degree more then - 0 - ",
            "Number of Rational Bezier curves  - ",
            "Number of  Bezier curves with less continuity - ",
            "Number of  other curves - ",
            "Number of BSpline pcurves with degree more then - 0 - ",
            "Number of BSpline pcurves with number of spans more then 0 - ",
            "Number of Rational BSpline pcurves - ",
            "Number of  BSpline pcurves with less continuity - ",
            "Number of Bezier pcurves with degree more then 0 -  - ",
            "Number of Rational Bezier pcurves  - ",
            "Number of  Bezier pcurves with less continuity - ",
            "Number of  other pcurves - ",
        ];
        Ok(labels.iter().zip(n).map(|(l, k)| format!("{l}{k}\n")).collect())
    }
    /// DRAW command `statshape`: counts of a shape's sub-shapes and kinds of geometry, in
    /// OCCT's report (the counts kentron gives; lines of none left out, as OCCT does).
    ///
    /// Example: `statshape result`
    pub fn statshape(&mut self, args: &[&str]) -> DrawResult {
        // `statshape shape`. "Curve Only C0" and "Surface Only C0" count edges and faces not
        // C1 inside them (by their curves' and surfaces' representation).
        need(args, 1, "statshape shape")?;
        use kentron::geometry::surface::Surface;
        let shape = self.shape_var(args[0])?;
        let m = &self.model;
        let unique = |s: Shape, k: NodeKind| m.unique_subshapes(s, k);
        let faces = unique(shape, NodeKind::Face);
        let edges = unique(shape, NodeKind::Edge);
        let in_faces: std::collections::HashSet<_> = faces.iter().flat_map(|&f| unique(f, NodeKind::Edge)).map(|e| e.node).collect();
        let in_shells: std::collections::HashSet<_> = unique(shape, NodeKind::Shell).iter().flat_map(|&s| unique(s, NodeKind::Face)).map(|f| f.node).collect();
        let wires = unique(shape, NodeKind::Wire);
        let wires_in_faces: std::collections::HashSet<_> = faces.iter().flat_map(|&f| unique(f, NodeKind::Wire)).map(|w| w.node).collect();
        // Edges as each face uses them (an edge used both ways counted twice).
        let mut oriented = std::collections::HashSet::new();
        let mut no_pcurve = 0;
        let mut trimmed2 = 0;
        for &face in &faces {
            for c in unique(face, NodeKind::CoEdge) {
                let id = kentron::topology::CoEdgeId::try_from(c.node).expect("a coedge");
                let edge = m[id].edge;
                oriented.insert((edge, m[id].orientation.compose(c.orientation) == kentron::topology::Orientation::Reversed));
                match m.effective_pcurve(id) {
                    None => no_pcurve += 1,
                    Some((Curve2::Trimmed(_), _, _)) => trimmed2 += 1,
                    _ => {}
                }
            }
        }
        let (mut bezier_surfaces, mut bspline_surfaces, mut c0_surfaces, mut offset_surfaces, mut trimmed_surfaces, mut big) = (0, 0, 0, 0, 0, 0);
        for &face in &faces {
            let Some(surface) = m.surface(face) else { continue };
            let (u0, u1, v0, v1) = m.face_uv_bounds(face);
            let (us, vs) = kentron::heal::surface_breaks(surface, u0, u1, v0, v1);
            c0_surfaces += (!us.is_empty() || !vs.is_empty()) as usize;
            match surface {
                Surface::Bezier(_) => bezier_surfaces += 1,
                Surface::BSpline(b) => {
                    bspline_surfaces += 1;
                    big += (b.poles().iter().map(|r| r.len()).sum::<usize>() > 8192) as usize;
                }
                Surface::Offset(_) => offset_surfaces += 1,
                Surface::Trimmed(_) => trimmed_surfaces += 1,
                _ => {}
            }
        }
        let (mut c0_curves, mut offset_curves, mut trimmed3) = (0, 0, 0);
        for &edge in &edges {
            let Some((curve, first, last)) = m.edge_curve(edge) else { continue };
            c0_curves += (!kentron::heal::curve_breaks(&curve, first, last).is_empty()) as usize;
            match curve {
                Curve::Offset(_) => offset_curves += 1,
                Curve::Trimmed(_) => trimmed3 += 1,
                Curve::BSpline(b) => big += (b.poles().len() > 8192) as usize,
                _ => {}
            }
        }
        let rows: Vec<(usize, &str)> = vec![
            (oriented.len(), "\t EDGE (Oriented)"),
            (edges.len(), "\t EDGE (Shared)"),
            (edges.iter().filter(|e| !in_faces.contains(&e.node)).count(), "\t EDGE (Free)"),
            (faces.len(), "\t FACE"),
            (faces.iter().filter(|f| !in_shells.contains(&f.node)).count(), "\t FACE (Free)"),
            (wires.iter().filter(|w| !wires_in_faces.contains(&w.node)).count(), "\t WIRE (Free)"),
            (unique(shape, NodeKind::Shell).len(), "\t SHELL"),
            (unique(shape, NodeKind::Solid).len(), "\t SOLID"),
            (unique(shape, NodeKind::Vertex).len(), "\t VERTEX (Shared)"),
            (wires.len(), "\t WIRE"),
            (faces.iter().filter(|&&f| unique(f, NodeKind::Wire).len() > 1).count(), "\tFace with more than one wire"),
            (no_pcurve, "\tNo pcurve"),
            (big, "\tbigspl : BSpline > 8192 poles"),
            (bezier_surfaces, "\tbezsur : BezierSurface"),
            (bspline_surfaces, "\tbspsur : BSplineSurface"),
            (c0_curves, "\tc0curv : Curve Only C0"),
            (c0_surfaces, "\tc0surf : Surface Only C0"),
            (offset_curves, "\tofcur  : Offset Curve(s)"),
            (offset_surfaces, "\tofsur  : Offset Surface"),
            (trimmed2, "\ttrc2d  : Trimmed Curve2d"),
            (trimmed3, "\ttrc3d  : Trimmed Curve3d"),
            (trimmed_surfaces, "\ttrimsu : RectangularTrimmedSurface"),
        ];
        let mut out = String::from("Count    Item\n-----    ----\n");
        for (k, label) in rows {
            if k > 0 {
                out += &format!("{k}{label}\n");
            }
        }
        Ok(out)
    }
    /// DRAW command `DT_SplitAngle`: 1 uses in 1 converted scripts.
    ///
    /// Example: `DT_SplitAngle result sh 90`
    pub fn DT_SplitAngle(&mut self, args: &[&str]) -> DrawResult {
        // `DT_SplitAngle result shape [max-angle=95]` (degrees, at least 1). The result is a
        // copy.
        if !(2..=3).contains(&args.len()) {
            return self.not_implemented("DT_SplitAngle", args);
        }
        let shape = self.shape_var(args[1])?;
        let max_angle = match args.get(2) { Some(a) => self.number(a)?.max(1.0), None => 95.0 };
        let copy = kentron::modify::copy(&mut self.model, shape, false).unwrap_or(shape);
        let split = kentron::heal::split_by_angle(&mut self.model, copy, max_angle.to_radians())?;
        self.set_shape_var(args[0], split)
    }
    /// DRAW command `DT_SplitByNumber`: 26 uses in 26 converted scripts.
    ///
    /// Example: `DT_SplitByNumber result a_1 10`
    pub fn DT_SplitByNumber(&mut self, args: &[&str]) -> DrawResult {
        // `DT_SplitByNumber result face number [number2]`: about `number` pieces of equal area,
        // or `number` along u by `number2` along v. The result is made from a copy.
        need(args, 3, "DT_SplitByNumber result face number [number2]")?;
        let face = self.shape_var(args[1])?;
        let count = |s: &str| -> Result<usize, DrawError> {
            let n = self.number(s)?;
            if n >= 1.0 { Ok(n as usize) } else { Err(DrawError::Failed(format!("incorrect number of parts: {s}"))) }
        };
        let parts = match args.get(3) {
            Some(second) => kentron::divide::Parts::Grid(count(args[2])?, count(second)?),
            None => kentron::divide::Parts::Total(count(args[2])?),
        };
        let copy = kentron::modify::copy(&mut self.model, face, false).unwrap_or(face);
        let divided = kentron::divide::divide_by_number(&mut self.model, copy, parts)?;
        self.set_shape_var(args[0], divided)
    }
    /// DRAW command `FAILED`: 1 uses in 1 converted scripts.
    ///
    /// Example: `FAILED /Not connected mesh inside face/ disconnected mesh`
    pub fn FAILED(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("FAILED", args)
    }
    /// DRAW command `GetShape`: 12 uses in 6 converted scripts.
    ///
    /// Example: `GetShape D 0:1:465:1:1:2 f`
    pub fn GetShape(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("GetShape", args)
    }
    /// DRAW command `OCC14376`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC14376 Face_8`
    pub fn OCC14376(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC14376", args)
    }
    /// DRAW command `OCC23429`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC23429 result a b`
    pub fn OCC23429(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC23429", args)
    }
    /// DRAW command `OCC23945`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC23945 s 0. -6.337359082e-009 x y z dux duy duz dvx dvy dvz`
    pub fn OCC23945(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC23945", args)
    }
    /// DRAW command `OCC23952intersect`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC23952intersect … a b`
    pub fn OCC23952intersect(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC23952intersect", args)
    }
    /// DRAW command `OCC23952sweep`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC23952sweep … a`
    pub fn OCC23952sweep(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC23952sweep", args)
    }
    /// DRAW command `OCC24667`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC24667 r w v3`
    pub fn OCC24667(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC24667", args)
    }
    /// DRAW command `OCC24834`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC24834`
    pub fn OCC24834(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC24834", args)
    }
    /// DRAW command `OCC26270`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26270 w2 result`
    pub fn OCC26270(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26270", args)
    }
    /// DRAW command `OCC26313`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26313 rr bb`
    pub fn OCC26313(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26313", args)
    }
    /// DRAW command `OCC26396`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26396 …`
    pub fn OCC26396(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26396", args)
    }
    /// DRAW command `OCC26446`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26446 r c1 c2`
    pub fn OCC26446(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26446", args)
    }
    /// DRAW command `OCC26485`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26485 a`
    pub fn OCC26485(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26485", args)
    }
    /// DRAW command `OCC26553`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC26553 …`
    pub fn OCC26553(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26553", args)
    }
    /// DRAW command `OCC26930`: 2 uses in 2 converted scripts.
    ///
    /// Example: `OCC26930 s c 0.0 14.3316`
    pub fn OCC26930(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC26930", args)
    }
    /// DRAW command `OCC28784`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC28784 result a`
    pub fn OCC28784(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC28784", args)
    }
    /// DRAW command `OCC30391`: 5 uses in 5 converted scripts.
    ///
    /// Example: `OCC30391 result a 0 0 20 20`
    pub fn OCC30391(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC30391", args)
    }
    /// DRAW command `OCC32744`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC32744 s`
    pub fn OCC32744(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC32744", args)
    }
    /// DRAW command `OCC363`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC363 OCC363Doc …`
    pub fn OCC363(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC363", args)
    }
    /// DRAW command `OCC367`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC367 result … … … … …`
    pub fn OCC367(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC367", args)
    }
    /// DRAW command `OCC369`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC369 result`
    pub fn OCC369(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC369", args)
    }
    /// DRAW command `OCC527`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC527 result`
    pub fn OCC527(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC527", args)
    }
    /// DRAW command `OCC566`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC566 result`
    pub fn OCC566(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC566", args)
    }
    /// DRAW command `OCC606`: 2 uses in 2 converted scripts.
    ///
    /// Example: `OCC606 result w2 -t`
    pub fn OCC606(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC606", args)
    }
    /// DRAW command `OCC7570`: 1 uses in 1 converted scripts.
    ///
    /// Example: `OCC7570 a`
    pub fn OCC7570(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("OCC7570", args)
    }
    /// DRAW command `Open`: 6 uses in 6 converted scripts.
    ///
    /// Example: `Open … D`
    pub fn Open(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("Open", args)
    }
    /// DRAW command `add`: 97 uses in 65 converted scripts.
    ///
    /// Example: `add r s`
    pub fn add(&mut self, args: &[&str]) -> DrawResult {
        // `add what where`: adds a shape to a container (in place).
        need(args, 2, "add what where")?;
        let (what, container) = (self.shape_var(args[0])?, self.shape_var(args[1])?);
        self.model.add(container, what)?;
        Ok(String::new())
    }
    /// DRAW command `addpolygonnode`: 1 uses in 1 converted scripts.
    ///
    /// Example: `addpolygonnode result 0 0`
    pub fn addpolygonnode(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("addpolygonnode", args)
    }
    /// DRAW command `addslide`: 9 uses in 9 converted scripts.
    ///
    /// Example: `addslide prism cont_1 E_3`
    pub fn addslide(&mut self, args: &[&str]) -> DrawResult {
        // `addslide kind edge face ...`: edges of the profile sliding on faces of the base,
        // for OCCT's gluing of the feature; kentron fuses or cuts the whole tool, so there is
        // nothing to record.
        if args.len() < 3 || args.len() % 2 == 0 {
            return self.not_implemented("addslide", args);
        }
        Ok(String::new())
    }
    /// DRAW command `addsweep`: 169 uses in 108 converted scripts.
    ///
    /// Example: `addsweep w`
    pub fn addsweep(&mut self, args: &[&str]) -> DrawResult {
        // `addsweep profile [vertex] [-T] [-R]`: a section (moved onto the spine with -T,
        // turned square to it with -R). A vertex of the spine it stands at is taken where
        // the profile lies there already (it stands where it lies, as without one).
        let options: Vec<&str> = args.iter().copied().filter(|a| a.starts_with('-')).collect();
        let names: Vec<&str> = args.iter().copied().filter(|a| !a.starts_with('-')).collect();
        if names.is_empty() || names.len() > 2 || options.iter().any(|o| *o != "-T" && *o != "-R") {
            return self.not_implemented("addsweep", args);
        }
        if let Some(vertex) = names.get(1) {
            let (v, profile) = (self.shape_var(vertex)?, self.shape_var(names[0])?);
            if v.kind() != NodeKind::Vertex {
                return self.not_implemented("addsweep", args);
            }
            let at = self.model.point(v);
            let gap = self.model.unique_subshapes(profile, NodeKind::Vertex).iter().map(|&w| self.model.point(w).distance(at)).fold(f64::INFINITY, f64::min);
            let size = properties::bounding_box(&self.model, profile);
            if !(gap <= (size.max - size.min).length()) {
                return self.not_implemented("addsweep", args);
            }
        }
        let profile = self.shape_var(names[0])?;
        self.sweep_profiles.push((profile, options.contains(&"-T"), options.contains(&"-R")));
        Ok(String::new())
    }
    /// DRAW command `appro`: 3 uses in 1 converted scripts.
    ///
    /// Example: `appro c1 16 ee1`
    pub fn appro(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("appro", args)
    }
    /// DRAW command `approxcurveonsurf`: 1 uses in 1 converted scripts.
    ///
    /// Example: `approxcurveonsurf res c2d sur`
    pub fn approxcurveonsurf(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("approxcurveonsurf", args)
    }
    /// DRAW command `appsurf`: 3 uses in 1 converted scripts.
    ///
    /// Example: `appsurf s1 c1 c2 c3`
    pub fn appsurf(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("appsurf", args)
    }
    /// DRAW command `baddcompound`: 6 uses in 5 converted scripts.
    ///
    /// Example: `baddcompound b1`
    pub fn baddcompound(&mut self, args: &[&str]) -> DrawResult {
        // `baddcompound c`: the compound's parts as objects for the general fuse.
        need(args, 1, "baddcompound compound")?;
        let c = self.shape_var(args[0])?;
        let parts: Vec<Shape> = self.model.children(c).collect();
        self.fuse_objects.extend(parts);
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `baddctools`: 4 uses in 4 converted scripts.
    ///
    /// Example: `baddctools c_2`
    pub fn baddctools(&mut self, args: &[&str]) -> DrawResult {
        // `baddctools c`: the compound's parts as tools for the general fuse.
        need(args, 1, "baddctools compound")?;
        let c = self.shape_var(args[0])?;
        let parts: Vec<Shape> = self.model.children(c).collect();
        self.fuse_tools.extend(parts);
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `bopcurves`: the curves where two faces meet.
    ///
    /// Example: `bopcurves f1 f2 -2d`
    pub fn bopcurves(&mut self, args: &[&str]) -> DrawResult {
        // `bopcurves f1 f2 [-2d|-2d1|-2d2] [-p u1 v1 u2 v2] [-v]`: the faces' section, a
        // curve per connected branch named c_1, c_2, ..., reported as OCCT does with the
        // farthest the curves stray from the faces' surfaces (at least 1e-7, OCCT's least).
        // (The pcurves and start points asked for are not needed to find them.)
        need(args, 2, "bopcurves f1 f2 [-2d|-2d1|-2d2] [-p u1 v1 u2 v2] [-v]")?;
        let (a, b) = (self.shape_var(args[0])?, self.shape_var(args[1])?);
        if a.kind() != NodeKind::Face || b.kind() != NodeKind::Face {
            return Err(DrawError::Failed(" Type mismatch F1".into()));
        }
        let section = boolean::section(&mut self.model, a, b)?;
        let edges: Vec<Shape> = self.model.unique_subshapes(section, NodeKind::Edge);
        let points: Vec<Shape> = self.model.unique_subshapes(section, NodeKind::Vertex).into_iter().filter(|v| !edges.iter().any(|&e| {
            let (p, q) = self.model.edge_vertices(e);
            p.node == v.node || q.node == v.node
        })).collect();
        if edges.is_empty() {
            return Ok(if points.is_empty() { " has no 3d curves\n has no 3d points\n".into() } else { format!(" has no 3d curves\n{} point(s) found.\n", points.len()) });
        }
        // The edges in connected chains, each in order.
        let chains = chain_edges(&self.model, &edges);
        let surfaces: Vec<Surface> = [a, b].iter().filter_map(|&f| self.model.surface(f).map(|s| s.transformed(&f.placement))).collect();
        let mut tolerance: f64 = 1e-7;
        let mut names = Vec::new();
        for (k, chain) in chains.iter().enumerate() {
            let mut samples = Vec::new();
            for &e in chain {
                let Some((c, first, last)) = self.model.edge_curve(e) else { continue };
                let (from, to) = if e.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
                let n = 32;
                for i in 0..=n {
                    if i == 0 && !samples.is_empty() {
                        continue;
                    }
                    samples.push(c.point(from + (to - from) * i as f64 / n as f64));
                }
            }
            for p in &samples {
                for s in &surfaces {
                    let (u, v, _) = kentron::geometry::extrema::nearest_on_surface(*p, s, None);
                    tolerance = tolerance.max(kentron::geometry::ParametricSurface::point(s, u, v).distance(*p));
                }
            }
            let curve = chain_curve(&self.model, chain)?;
            let name = format!("c_{}", k + 1);
            self.set(&name, DrawObject::Curve(curve));
            names.push(name);
        }
        Ok(format!("Tolerance Reached={}\n{} curve(s) found.\n{}\n", format_real(tolerance), chains.len(), names.join(" ")))
    }
    /// DRAW command `baddobjects`: 269 uses in 255 converted scripts.
    ///
    /// Example: `baddobjects a b`
    pub fn baddobjects(&mut self, args: &[&str]) -> DrawResult {
        // `baddobjects s1 s2 ...`: objects for the general fuse.
        for name in args {
            let shape = self.shape_var(name)?;
            self.fuse_objects.push(shape);
        }
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `baddtools`: 142 uses in 133 converted scripts.
    ///
    /// Example: `baddtools f3`
    pub fn baddtools(&mut self, args: &[&str]) -> DrawResult {
        // `baddtools s1 s2 ...`: tools for the general fuse.
        for name in args {
            let shape = self.shape_var(name)?;
            self.fuse_tools.push(shape);
        }
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `bapibop`: 10 uses in 5 converted scripts.
    ///
    /// Example: `bapibop result 4`
    pub fn bapibop(&mut self, args: &[&str]) -> DrawResult {
        // `bapibop result operation`: `bbop` through OCCT's API classes, the result then
        // simplified as `bsimplify` asks.
        need(args, 2, "bapibop result operation")?;
        self.bbop(args)?;
        let (edges, faces) = self.simplify;
        if edges || faces {
            let result = self.shape_var(args[0])?;
            let options = kentron::unify::UnifyOptions { edges, faces, ..Default::default() };
            let simplified = kentron::unify::unify_same_domain_with(&mut self.model, result, options)?;
            return self.set_shape_var(args[0], simplified);
        }
        Ok(String::new())
    }
    /// DRAW command `bapibuild`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bapibuild result`
    pub fn bapibuild(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("bapibuild", args)
    }
    /// DRAW command `bapisplit`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bapisplit result`
    pub fn bapisplit(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("bapisplit", args)
    }
    /// DRAW command `bbop`: 95 uses in 72 converted scripts.
    ///
    /// Example: `bbop result 4`
    pub fn bbop(&mut self, args: &[&str]) -> DrawResult {
        // `bbop result operation`: 0 common, 1 fuse, 2 cut, 3 the tools cut by the objects,
        // 4 section, of the objects and the tools, every argument split by the others.
        need(args, 2, "bbop result operation")?;
        // (By number, or by name.)
        let operation: i32 = match args[1].to_ascii_lowercase().as_str() {
            "common" => 0,
            "fuse" => 1,
            "cut" => 2,
            "tuc" | "cut21" => 3,
            "section" => 4,
            n => n.parse().map_err(|_| DrawError::Failed(format!("bbop: bad operation {}", args[1])))?,
        };
        use boolean::Combination as C;
        let combination = match operation {
            0 => C::Intersection,
            1 => C::Union,
            2 => C::Difference,
            3 => C::ReversedDifference,
            4 => C::Section,
            _ => return Err(DrawError::Failed(format!("bbop: bad operation {operation}"))),
        };
        if self.fuse_objects.is_empty() {
            return Err(DrawError::Failed("bbop: no arguments".into()));
        }
        let (objects, tools) = (self.fuse_objects.clone(), self.fuse_tools.clone());
        let result = boolean::combine(&mut self.model, &objects, &tools, combination)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bbuild`: 110 uses in 104 converted scripts.
    ///
    /// Example: `bbuild r`
    pub fn bbuild(&mut self, args: &[&str]) -> DrawResult {
        // `bbuild result [-t]`: the general fuse: every argument split by the others, the
        // pieces together.
        need(args, 1, "bbuild result")?;
        let start = std::time::Instant::now();
        let all: Vec<Shape> = self.fuse_cells()?.iter().map(|c| c.shape).collect();
        let result = self.compound_of(&all)?;
        self.set_shape_var(args[0], result)?;
        // (With `-t`, the time taken, as OCCT reports it.)
        Ok(if args.contains(&"-t") { format!(" Tps: {}\n", start.elapsed().as_secs_f64()) } else { String::new() })
    }
    /// DRAW command `bcadd`: 120 uses in 72 converted scripts.
    ///
    /// Example: `bcadd result b 1 -m 1 -u`
    pub fn bcadd(&mut self, args: &[&str]) -> DrawResult {
        // `bcadd result s1 1 s2 0 ... [-m material [-u]]`: the cells inside the shapes
        // marked 1 and outside those marked 0 added, with the material (merged with its
        // other cells, -u).
        let (selected, material, update) = self.cells_selection("bcadd", args)?;
        for i in selected {
            if !self.cells_taken.iter().any(|&(j, _)| j == i) {
                self.cells_taken.push((i, material));
            }
        }
        if update && material != 0 {
            self.cells_merged.insert(material);
        }
        self.set_cells_result(args[0])
    }
    /// DRAW command `bcaddall`: 13 uses in 13 converted scripts.
    ///
    /// Example: `bcaddall result`
    pub fn bcaddall(&mut self, args: &[&str]) -> DrawResult {
        // `bcaddall result [-m material [-u]]`: every cell taken.
        need(args, 1, "bcaddall result [-m material [-u]]")?;
        let mut material = 0;
        if let Some(k) = args.iter().position(|a| *a == "-m") {
            material = args.get(k + 1).and_then(|m| m.parse().ok()).ok_or_else(|| DrawError::Failed("bcaddall: bad material".into()))?;
        }
        let count = self.fuse_cells()?.len();
        self.cells_taken = (0..count).map(|i| (i, material)).collect();
        self.cells_merged.clear();
        if material != 0 && args.contains(&"-u") {
            self.cells_merged.insert(material);
        }
        self.set_cells_result(args[0])
    }
    /// DRAW command `bcbuild`: 78 uses in 78 converted scripts.
    ///
    /// Example: `bcbuild rx`
    pub fn bcbuild(&mut self, args: &[&str]) -> DrawResult {
        // `bcbuild result`: the cells of the general fuse's arguments, all taken.
        need(args, 1, "bcbuild result")?;
        let count = self.fuse_cells()?.len();
        self.cells_taken = (0..count).map(|i| (i, 0)).collect();
        self.cells_merged.clear();
        self.set_cells_result(args[0])
    }
    /// DRAW command `bclassify`: 2 uses in 2 converted scripts.
    ///
    /// Example: `bclassify b1 p`
    pub fn bclassify(&mut self, args: &[&str]) -> DrawResult {
        // `bclassify solid point [tolerance]`: where the point lies.
        need(args, 2, "bclassify solid point [tolerance]")?;
        let solid = self.shape_var(args[0])?;
        let Some(DrawObject::Point(p)) = self.get(args[1]).cloned() else {
            return Err(DrawError::Failed(format!("{} is not a point", args[1])));
        };
        let place = match boolean::classify_point(&self.model, solid, p)? {
            boolean::PointLocation::Inside => "IN",
            boolean::PointLocation::Outside => "OUT",
            _ => "ON",
        };
        Ok(format!("The point is {place} shape"))
    }
    /// DRAW command `bclear`: 2 uses in 2 converted scripts.
    ///
    /// Example: `bclear`
    pub fn bclear(&mut self, _args: &[&str]) -> DrawResult {
        // `bclear`: no arguments for the general fuse.
        self.fuse_objects.clear();
        self.fuse_tools.clear();
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `bclearobjects`: 269 uses in 255 converted scripts.
    ///
    /// Example: `bclearobjects`
    pub fn bclearobjects(&mut self, _args: &[&str]) -> DrawResult {
        // `bclearobjects`: no objects for the general fuse.
        self.fuse_objects.clear();
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `bcleartools`: 263 uses in 249 converted scripts.
    ///
    /// Example: `bcleartools`
    pub fn bcleartools(&mut self, _args: &[&str]) -> DrawResult {
        // `bcleartools`: no tools for the general fuse.
        self.fuse_tools.clear();
        self.fuse_cells = None;
        Ok(String::new())
    }
    /// DRAW command `bcmakecontainers`: 10 uses in 10 converted scripts.
    ///
    /// Example: `bcmakecontainers result`
    pub fn bcmakecontainers(&mut self, args: &[&str]) -> DrawResult {
        // `bcmakecontainers result`: the cells as they are (OCCT groups touching ones into
        // containers, which changes no measure).
        need(args, 1, "bcmakecontainers result")?;
        self.set_cells_result(args[0])
    }
    /// DRAW command `bcommon`: 299 uses in 294 converted scripts.
    ///
    /// Example: `bcommon result b a`
    pub fn bcommon(&mut self, args: &[&str]) -> DrawResult {
        // `bcommon result object tool`.
        if args.len() != 3 {
            return self.not_implemented("bcommon", args);
        }
        let (object, tool) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let result = boolean::intersection(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bcremove`: 14 uses in 10 converted scripts.
    ///
    /// Example: `bcremove result s1 1 s2 1 s3 1`
    pub fn bcremove(&mut self, args: &[&str]) -> DrawResult {
        // `bcremove result s1 1 s2 0 ...`: the cells inside the shapes marked 1 and outside
        // those marked 0 taken out.
        let (selected, _, _) = self.cells_selection("bcremove", args)?;
        self.cells_taken.retain(|(i, _)| !selected.contains(i));
        self.set_cells_result(args[0])
    }
    /// DRAW command `bcremoveall`: 74 uses in 59 converted scripts.
    ///
    /// Example: `bcremoveall`
    pub fn bcremoveall(&mut self, _args: &[&str]) -> DrawResult {
        // `bcremoveall`: no cells taken.
        self.cells_taken.clear();
        self.cells_merged.clear();
        Ok(String::new())
    }
    /// DRAW command `bcremoveint`: 10 uses in 10 converted scripts.
    ///
    /// Example: `bcremoveint result`
    pub fn bcremoveint(&mut self, args: &[&str]) -> DrawResult {
        // `bcremoveint result`: the cells of each material merged, their boundaries between
        // them gone.
        need(args, 1, "bcremoveint result")?;
        let materials: Vec<i64> = self.cells_taken.iter().map(|&(_, m)| m).filter(|&m| m != 0).collect();
        self.cells_merged.extend(materials);
        self.set_cells_result(args[0])
    }
    /// DRAW command `bcut`: 735 uses in 583 converted scripts.
    ///
    /// Example: `bcut s s1 s2`
    pub fn bcut(&mut self, args: &[&str]) -> DrawResult {
        // `bcut result object tool`.
        if args.len() != 3 {
            return self.not_implemented("bcut", args);
        }
        let (object, tool) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let result = boolean::difference(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bcutblend`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bcutblend result box1 box2 0.01`
    pub fn bcutblend(&mut self, args: &[&str]) -> DrawResult {
        // `bcutblend result shape1 shape2 radius [-d]`: the boolean operation, then its section
        // edges filleted (the boolean's result where that fails, like OCCT).
        need(args, 4, "bcutblend result shape1 shape2 radius [-d]")?;
        let (a, b) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let radius = self.number(args[3])?;
        let (shape, sections) = boolean::difference_with_sections(&mut self.model, a, b)?;
        if args.get(4) == Some(&"-d") {
            self.set_shape_var(&format!("{}_bop", args[0]), shape)?;
        }
        let edges: Vec<(Shape, f64)> = sections.into_iter().map(|e| (e, radius)).collect();
        let result = if edges.is_empty() { shape } else { blend::fillet_edges(&mut self.model, shape, &edges).unwrap_or(shape) };
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `beziercurve`: 88 uses in 80 converted scripts.
    ///
    /// Example: `beziercurve w1 5 0 0 0 20 0 0 20 5 0 25 10 0 10 20 0`
    pub fn beziercurve(&mut self, args: &[&str]) -> DrawResult {
        // `beziercurve name n x1 y1 z1 [w1] ... xn yn zn [wn]`.
        need(args, 2, "beziercurve name n poles...")?;
        let n = self.number(args[1])? as usize;
        let v = self.numbers(&args[2..])?;
        let stride = if v.len() == 4 * n { 4 } else if v.len() == 3 * n { 3 } else {
            return Err(DrawError::Failed("beziercurve: wrong number of coordinates".into()));
        };
        let poles: Vec<Point> = v.chunks(stride).map(|c| Point::new(c[0], c[1], c[2])).collect();
        let curve = if stride == 4 {
            BezierCurve::rational(poles, v.chunks(4).map(|c| c[3]).collect())
        } else {
            BezierCurve::new(poles)
        };
        self.set(args[0], DrawObject::Curve(Curve::Bezier(curve)));
        Ok(String::new())
    }
    /// DRAW command `beziersurf`: 17 uses in 17 converted scripts.
    ///
    /// Example: `beziersurf bz 3 3 0 0 0 1 -1 0 2 0 0 0 -1 1 1 0 1 2 -1 1 0 0 2 1 -1 2 2 0 2`
    pub fn beziersurf(&mut self, args: &[&str]) -> DrawResult {
        // `beziersurf name nu nv poles (x y z [w]), u varying fastest`.
        need(args, 3, "beziersurf name nu nv poles...")?;
        let (nu, nv) = (self.number(args[1])? as usize, self.number(args[2])? as usize);
        let v = self.numbers(&args[3..])?;
        let stride = match nu * nv {
            0 => 0,
            n if v.len() == 4 * n => 4,
            n if v.len() == 3 * n => 3,
            _ => 0,
        };
        if stride == 0 {
            return Err(DrawError::Failed("beziersurf: wrong number of coordinates".into()));
        }
        let at = |i: usize, j: usize| &v[(j * nu + i) * stride..][..stride];
        let poles: Vec<Vec<Point>> = (0..nu).map(|i| (0..nv).map(|j| { let c = at(i, j); Point::new(c[0], c[1], c[2]) }).collect()).collect();
        let surface = BezierSurface::try_new(poles).ok_or_else(|| DrawError::Failed("beziersurf: invalid pole grid".into()))?;
        let surface = match stride {
            4 => BezierSurface::rational(surface.poles().to_vec(), (0..nu).map(|i| (0..nv).map(|j| at(i, j)[3]).collect()).collect()),
            _ => surface,
        };
        self.set(args[0], DrawObject::Surface(Surface::Bezier(surface)));
        Ok(String::new())
    }
    /// DRAW command `bfillds`: 269 uses in 255 converted scripts.
    ///
    /// Example: `bfillds`
    pub fn bfillds(&mut self, _args: &[&str]) -> DrawResult {
        // `bfillds [-t]`: the arguments are intersected when a result is asked for.
        if self.fuse_objects.is_empty() && self.fuse_tools.is_empty() {
            return Ok("No objects to process\n".into());
        }
        Ok(String::new())
    }
    /// DRAW command `bfuse`: 630 uses in 538 converted scripts.
    ///
    /// Example: `bfuse s c b`
    pub fn bfuse(&mut self, args: &[&str]) -> DrawResult {
        // `bfuse result object tool`.
        if args.len() != 3 {
            return self.not_implemented("bfuse", args);
        }
        let (object, tool) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let result = boolean::union(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bfuseblend`: 18 uses in 18 converted scripts.
    ///
    /// Example: `bfuseblend result s b 1*SCALE1`
    pub fn bfuseblend(&mut self, args: &[&str]) -> DrawResult {
        // `bfuseblend result shape1 shape2 radius [-d]`: the boolean operation, then its section
        // edges filleted (the boolean's result where that fails, like OCCT).
        need(args, 4, "bfuseblend result shape1 shape2 radius [-d]")?;
        let (a, b) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let radius = self.number(args[3])?;
        let (shape, sections) = boolean::union_with_sections(&mut self.model, a, b)?;
        if args.get(4) == Some(&"-d") {
            self.set_shape_var(&format!("{}_bop", args[0]), shape)?;
        }
        let edges: Vec<(Shape, f64)> = sections.into_iter().map(|e| (e, radius)).collect();
        let result = if edges.is_empty() { shape } else { blend::fillet_edges(&mut self.model, shape, &edges).unwrap_or(shape) };
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bfuzzyvalue`: 26 uses in 21 converted scripts.
    ///
    /// Example: `bfuzzyvalue 5.e-6`
    pub fn bfuzzyvalue(&mut self, args: &[&str]) -> DrawResult {
        // `bfuzzyvalue value`: OCCT's booleans then take gaps up to the value as touching.
        // Kentron's have no such setting (they meet within their own tolerances): accepted,
        // the script's checks telling whether the results agree.
        need(args, 1, "bfuzzyvalue value")?;
        args[0].parse::<f64>().map_err(|_| DrawError::Failed(format!("bfuzzyvalue: bad value {}", args[0])))?;
        Ok(String::new())
    }
    /// DRAW command `bglue`: 2 uses in 2 converted scripts.
    ///
    /// Example: `bglue 1`
    pub fn bglue(&mut self, args: &[&str]) -> DrawResult {
        // `bglue 0|1|2`: tells OCCT's booleans that the arguments only touch (sharing faces),
        // for speed; kentron's results are the same without it.
        need(args, 1, "bglue 0|1|2")?;
        Ok(String::new())
    }
    /// DRAW command `bhaspc`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bhaspc ac a do`
    pub fn bhaspc(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("bhaspc", args)
    }
    /// DRAW command `binrestore`: 25 uses in 17 converted scripts.
    ///
    /// Example: `binrestore … a`
    pub fn binrestore(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("binrestore", args)
    }
    /// DRAW command `blend`: 354 uses in 337 converted scripts.
    ///
    /// Example: `blend result s 10 s_5`
    pub fn blend(&mut self, args: &[&str]) -> DrawResult {
        // `blend result shape radius edge [radius edge ...]`: fillets (the trailing shape
        // options of OCCT's command are not supported).
        if args.len() < 4 || args.len() % 2 != 0 {
            return self.not_implemented("blend", args);
        }
        let shape = self.shape_var(args[1])?;
        let mut edges = Vec::new();
        for pair in args[2..].chunks(2) {
            edges.push((self.shape_var(pair[1])?, self.number(pair[0])?));
        }
        let result = blend::fillet_edges(&mut self.model, shape, &edges)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bnondestructive`: 4 uses in 4 converted scripts.
    ///
    /// Example: `bnondestructive 1`
    pub fn bnondestructive(&mut self, args: &[&str]) -> DrawResult {
        // `bnondestructive 0|1`: whether booleans may modify their operands; kentron's
        // never do.
        need(args, 1, "bnondestructive 0|1")?;
        Ok(String::new())
    }
    /// DRAW command `bop`: 2746 uses in 2685 converted scripts.
    ///
    /// Example: `bop a b`
    pub fn bop(&mut self, args: &[&str]) -> DrawResult {
        // `bop object tool`: prepares a boolean operation (run by `bopfuse`, `bopcut`, ...).
        need(args, 2, "bop object tool")?;
        self.bop_operands = Some((self.shape_var(args[0])?, self.shape_var(args[1])?));
        Ok(String::new())
    }
    /// DRAW command `bopaddpcs`: 4 uses in 2 converted scripts.
    ///
    /// Example: `bopaddpcs a`
    pub fn bopaddpcs(&mut self, args: &[&str]) -> DrawResult {
        // `bopaddpcs shape`: pcurves for the edges of its planar faces that lack them.
        need(args, 1, "bopaddpcs shape")?;
        let shape = self.shape_var(args[0])?;
        kentron::heal::add_plane_pcurves(&mut self.model, shape)?;
        Ok(String::new())
    }
    /// DRAW command `bopargcheck`: 20 uses in 16 converted scripts.
    ///
    /// Example: `bopargcheck b_1 b_2`
    pub fn bopargcheck(&mut self, args: &[&str]) -> DrawResult {
        // `bopargcheck shape1 [shape2] [-F|O|C|T|S|U] [/options] [#BF]`: whether the shapes
        // can be arguments of a boolean: valid (`checkshape`) and without sub-shapes meeting
        // outside their shared vertices (`bopcheck`). The operation and test options are
        // not taken into account.
        need(args, 1, "bopargcheck shape1 [shape2] [options]")?;
        let shapes: Vec<Shape> = args.iter().filter(|a| !a.starts_with(['-', '/', '#'])).map(|a| self.shape_var(a)).collect::<DrawResult<_>>()?;
        let mut faulty = false;
        for &shape in &shapes {
            faulty |= !kentron::validity::check(&self.model, shape).is_valid();
            faulty |= !kentron::validity::self_interferences(&mut self.model, shape)?.is_empty();
        }
        Ok(if faulty {
            "Faulties, that can not be treated by BOP, are detected.\n".into()
        } else {
            "Shape(s) seem(s) to be valid for BOP.\n".into()
        })
    }
    /// DRAW command `bopbsolid`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bopbsolid r s`
    pub fn bopbsolid(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("bopbsolid", args)
    }
    /// DRAW command `bopcheck`: 34 uses in 25 converted scripts.
    ///
    /// Example: `bopcheck b_1`
    pub fn bopcheck(&mut self, args: &[&str]) -> DrawResult {
        // `bopcheck shape [level] [-t]`: the pairs of its sub-shapes that meet without sharing
        // a vertex (`kentron::validity::self_interferences`), up to the level's kinds (0 V/V,
        // 1 V/E, 2 E/E, 3 V/F, 4 E/F, 5 and above F/F; kentron does not check solids), each
        // as `x<n> x<m>` with the shapes so named (numbered as OCCT's checker numbers them:
        // the shape, then its sub-shapes depth first).
        need(args, 1, "bopcheck shape [level] [-t]")?;
        let shape = self.shape_var(args[0])?;
        let level: usize = match args.get(1).filter(|a| !a.starts_with('-')) {
            Some(l) => l.parse().map_err(|_| DrawError::Failed(format!("bopcheck: bad level {l}")))?,
            None => 9,
        };
        let numbers: std::collections::HashMap<NumberingKey, usize> = occt_numbering(&self.model, shape).into_iter().enumerate().map(|(i, s)| (numbering_key(s), i)).collect();
        let number = |s: Shape| numbers.get(&numbering_key(s)).copied().unwrap_or(usize::MAX);
        let rank = |k: NodeKind| match k {
            NodeKind::Vertex => 0,
            NodeKind::Edge => 1,
            _ => 2,
        };
        const KINDS: [&str; 6] = ["V/V", "V/E", "E/E", "V/F", "E/F", "F/F"];
        let mut found: Vec<(usize, usize, usize, Shape, Shape)> = Vec::new();
        for i in kentron::validity::self_interferences(&mut self.model, shape)? {
            let (a, b) = (rank(i.first.kind()), rank(i.second.kind()));
            let kind = match (a, b) {
                (0, 0) => 0,
                (0, 1) => 1,
                (1, 1) => 2,
                (0, 2) => 3,
                (1, 2) => 4,
                _ => 5,
            };
            if kind <= level {
                found.push((kind, number(i.first), number(i.second), i.first, i.second));
            }
        }
        found.sort_by_key(|f| (f.0, f.1, f.2));
        let mut out = String::new();
        for (kind, n1, n2, first, second) in found.iter().copied() {
            self.set_shape_var(&format!("x{n1}"), first)?;
            self.set_shape_var(&format!("x{n2}"), second)?;
            out.push_str(&format!("{}: x{n1} x{n2} \n", KINDS[kind]));
        }
        if found.is_empty() {
            out.push_str(" This shape seems to be OK.\n");
        }
        Ok(out)
    }
    /// DRAW command `bopcommon`: 709 uses in 706 converted scripts.
    ///
    /// Example: `bopcommon result`
    pub fn bopcommon(&mut self, args: &[&str]) -> DrawResult {
        // `bopcommon result`: the common of the last `bop`.
        need(args, 1, "bopcommon result")?;
        let (object, tool) = self.bop_operands.ok_or_else(|| DrawError::Failed("bop has not been run".into()))?;
        let result = boolean::intersection(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bopcut`: 630 uses in 617 converted scripts.
    ///
    /// Example: `bopcut result`
    pub fn bopcut(&mut self, args: &[&str]) -> DrawResult {
        // `bopcut result`: the cut of the last `bop`.
        need(args, 1, "bopcut result")?;
        let (object, tool) = self.bop_operands.ok_or_else(|| DrawError::Failed("bop has not been run".into()))?;
        let result = boolean::difference(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bopfuse`: 635 uses in 600 converted scripts.
    ///
    /// Example: `bopfuse result`
    pub fn bopfuse(&mut self, args: &[&str]) -> DrawResult {
        // `bopfuse result`: the fuse of the last `bop`.
        need(args, 1, "bopfuse result")?;
        let (object, tool) = self.bop_operands.ok_or_else(|| DrawError::Failed("bop has not been run".into()))?;
        let result = boolean::union(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `boporigin`: 4 uses in 1 converted scripts.
    ///
    /// Example: `boporigin r_1`
    pub fn boporigin(&mut self, args: &[&str]) -> DrawResult {
        // `boporigin shape`: the sub-shapes of the general fuse's arguments that `shape` (a
        // piece of its result) comes from, as `shape_or` (one, or a compound): those of its
        // kind holding its points (a vertex's point; an edge's ends and middle; a face's or a
        // solid's centre of mass, which must lie on it).
        need(args, 1, "boporigin shape")?;
        let shape = self.shape_var(args[0])?;
        let kind = shape.kind();
        let on = |draw: &mut Self, p: kentron::math::Point, s: Shape| -> DrawResult<bool> {
            let v = kentron::make::vertex(&mut draw.model, p);
            Ok(kentron::distance::distance(&mut draw.model, v, s)?.value <= 1e-6)
        };
        let points: Vec<kentron::math::Point> = match kind {
            NodeKind::Vertex => vec![self.model.point(shape)],
            NodeKind::Edge => {
                let (start, end) = self.model.edge_vertices(shape);
                let middle = kentron::properties::length(&self.model, shape).center_of_mass;
                vec![self.model.point(start), self.model.point(end), middle]
            }
            NodeKind::Face => vec![kentron::properties::area(&self.model, shape).center_of_mass],
            NodeKind::Solid => vec![kentron::properties::volume(&self.model, shape).center_of_mass],
            _ => return self.not_implemented("boporigin", args),
        };
        // The point must lie on the shape itself (the centre of a curved edge or face does not).
        let inside = |draw: &mut Self, p: kentron::math::Point, s: Shape| -> DrawResult<bool> {
            if s.kind() == NodeKind::Solid {
                return Ok(boolean::classify_point(&draw.model, s, p)? == boolean::PointLocation::Inside);
            }
            on(draw, p, s)
        };
        for &p in &points {
            if !inside(self, p, shape)? {
                return self.not_implemented("boporigin", args);
            }
        }
        let arguments: Vec<Shape> = self.fuse_objects.iter().chain(&self.fuse_tools).copied().collect();
        let mut origins = Vec::new();
        for argument in arguments {
            for candidate in self.model.unique_subshapes(argument, kind) {
                let mut holds = true;
                for &p in &points {
                    if !inside(self, p, candidate)? {
                        holds = false;
                        break;
                    }
                }
                if holds {
                    origins.push(candidate);
                }
            }
        }
        let name = format!("{}_or", args[0]);
        match origins.as_slice() {
            [] => Ok(" no origins found\n".into()),
            [one] => {
                self.set_shape_var(&name, *one)?;
                Ok(format!("1 origin found\n{name}\n"))
            }
            _ => {
                let compound = self.compound_of(&origins)?;
                self.set_shape_var(&name, compound)?;
                Ok(format!("{} origins found\n{name}\n", origins.len()))
            }
        }
    }
    /// DRAW command `bopsection`: 152 uses in 151 converted scripts.
    ///
    /// Example: `bopsection result`
    pub fn bopsection(&mut self, args: &[&str]) -> DrawResult {
        // `bopsection result`: the section of the last `bop`.
        need(args, 1, "bopsection result")?;
        let (object, tool) = self.bop_operands.ok_or_else(|| DrawError::Failed("bop has not been run".into()))?;
        let result = boolean::section(&mut self.model, object, tool)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `boptuc`: 622 uses in 622 converted scripts.
    ///
    /// Example: `boptuc result`
    pub fn boptuc(&mut self, args: &[&str]) -> DrawResult {
        // `boptuc result`: the tool minus the object of the last `bop`.
        need(args, 1, "boptuc result")?;
        let (object, tool) = self.bop_operands.ok_or_else(|| DrawError::Failed("bop has not been run".into()))?;
        let result = boolean::difference(&mut self.model, tool, object)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `bounding`: 19 uses in 15 converted scripts.
    ///
    /// Example: `bounding sh1 -dump -shape result`
    pub fn bounding(&mut self, args: &[&str]) -> DrawResult {
        // `bounding shape [-save xmin ymin zmin xmax ymax zmax] [-dump|-print] [-shape name]`
        // (`-optimal`, `-exttoler` and `-noTriangulation` accepted; kentron's boxes are tight).
        need(args, 1, "bounding shape [options]")?;
        // (Or the box given by its corners, `xmin ymin zmin xmax ymax zmax`.)
        let corners = args.len() >= 6 && args[..6].iter().all(|a| a.parse::<f64>().is_ok());
        let (b, mut k) = if corners {
            let v = self.numbers(&args[..6])?;
            (kentron::math::BoundingBox::new(Point::new(v[0], v[1], v[2]), Point::new(v[3], v[4], v[5])), 6)
        } else {
            let shape = self.shape_var(args[0])?;
            (properties::tight_bounding_box(&self.model, shape), 1)
        };
        let values = [b.min.x, b.min.y, b.min.z, b.max.x, b.max.y, b.max.z];
        while k < args.len() {
            match args[k] {
                "-save" => {
                    for (name, value) in args.get(k + 1..k + 7).unwrap_or(&[]).iter().zip(values) {
                        self.set(name, DrawObject::Number(value));
                    }
                    k += 7;
                }
                "-shape" => {
                    let name = args.get(k + 1).copied().ok_or_else(|| DrawError::Failed("bounding: -shape needs a name".into()))?;
                    let size = b.max - b.min;
                    let tiny = |x: f64| x.max(1e-7);
                    let solid = Cuboid::from_corners(b.min, b.min + Vector::new(tiny(size.x), tiny(size.y), tiny(size.z))).build(&mut self.model)?;
                    self.set_shape_var(name, solid)?;
                    k += 2;
                }
                "-dump" | "-print" | "-optimal" | "-exttoler" | "-noTriangulation" | "-obb" => k += 1,
                other => return Err(DrawError::Failed(format!("bounding: unknown option {other}"))),
            }
        }
        Ok(values.iter().map(|&v| format_real(v)).collect::<Vec<_>>().join(" "))
    }
    /// DRAW command `bounds`: 6 uses in 5 converted scripts.
    ///
    /// Example: `bounds c2d2 u1 u2`
    pub fn bounds(&mut self, args: &[&str]) -> DrawResult {
        // `bounds curve u1 u2` or `bounds surface u1 u2 v1 v2`: the parameter ranges, into
        // variables (infinite ones as ±2e100, as OCCT's).
        need(args, 3, "bounds curve|surface u1 u2 [v1 v2]")?;
        let finite = |x: f64| x.clamp(-2e100, 2e100);
        let values = match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(c)) => vec![c.first_parameter(), c.last_parameter()],
            Some(DrawObject::Curve2(c)) => vec![c.first_parameter(), c.last_parameter()],
            Some(DrawObject::Surface(s)) => {
                let (u0, u1, v0, v1) = kentron::geometry::ParametricSurface::bounds(&s);
                vec![u0, u1, v0, v1]
            }
            _ => return Err(DrawError::Failed(format!("bounds: {} is not a curve or surface", args[0]))),
        };
        for (name, value) in args[1..].iter().zip(values) {
            self.set(name, DrawObject::Number(finite(value)));
        }
        Ok(String::new())
    }
    /// DRAW command `box`: 2511 uses in 1659 converted scripts.
    ///
    /// Example: `box s 100 100 100`
    pub fn r#box(&mut self, args: &[&str]) -> DrawResult {
        // `box name [x y z] dx dy dz` (the `-min`, `-max`, `-size`, `-dir`, `-xdir` and
        // `-preview` forms are not supported).
        need(args, 4, "box name [x y z] dx dy dz")?;
        const OPTIONS: [&str; 6] = ["-min", "-max", "-size", "-dir", "-xdir", "-preview"];
        if args[1..].iter().any(|a| OPTIONS.contains(&a.to_ascii_lowercase().as_str())) {
            return self.not_implemented("box", args);
        }
        let values = args[1..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let (origin, size) = match values.len() {
            3 => (Point::ORIGIN, &values[..]),
            6 => (Point::new(values[0], values[1], values[2]), &values[3..]),
            _ => return Err(DrawError::Failed("box: syntax error".into())),
        };
        // As BRepPrimAPI_MakeBox, a negative size extends the box backwards from the corner.
        let far = origin + Vector::new(size[0], size[1], size[2]);
        if size.iter().any(|&d| d.abs() <= tolerance::LINEAR) {
            return Err(DrawError::Failed("box: null dimension".into()));
        }
        let shape = Cuboid::from_corners(origin, far).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `breducetolerance`: 4 uses in 3 converted scripts.
    ///
    /// Example: `breducetolerance a`
    pub fn breducetolerance(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("breducetolerance", args)
    }
    /// DRAW command `brestore`: 23 uses in 19 converted scripts.
    ///
    /// Example: `brestore … b1`
    pub fn brestore(&mut self, args: &[&str]) -> DrawResult {
        // `brestore file name`: the shape of a BRep file (as `restore` reads one).
        need(args, 2, "brestore file name")?;
        self.restore(args)
    }
    /// DRAW command `brollingball`: 4 uses in 4 converted scripts.
    ///
    /// Example: `brollingball result s 15 @ f1 f2`
    pub fn brollingball(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("brollingball", args)
    }
    /// DRAW command `brunparallel`: 2 uses in 2 converted scripts.
    ///
    /// Example: `brunparallel 1`
    pub fn brunparallel(&mut self, args: &[&str]) -> DrawResult {
        // `brunparallel 0|1`: whether OCCT's booleans run in parallel; kentron's choose for
        // themselves, with the same results.
        need(args, 1, "brunparallel 0|1")?;
        Ok(String::new())
    }
    /// DRAW command `bsection`: 219 uses in 211 converted scripts.
    ///
    /// Example: `bsection sec part tool`
    pub fn bsection(&mut self, args: &[&str]) -> DrawResult {
        // `bsection result shape1 shape2` (options approximating or computing pcurves do not
        // change the edges).
        let operands: Vec<&str> = args.iter().copied().filter(|a| !a.starts_with('-')).collect();
        if operands.len() != 3 {
            return self.not_implemented("bsection", args);
        }
        let (a, b) = (self.shape_var(operands[1])?, self.shape_var(operands[2])?);
        let result = boolean::section(&mut self.model, a, b)?;
        self.set_shape_var(operands[0], result)
    }
    /// DRAW command `bsimplify`: 5 uses in 3 converted scripts.
    ///
    /// Example: `bsimplify -e 1`
    pub fn bsimplify(&mut self, args: &[&str]) -> DrawResult {
        // `bsimplify [-e 0|1] [-f 0|1] [-a angle]`: whether `bapibop` unifies the result's
        // edges and faces lying on one curve or surface (the angle is not taken into account).
        let mut rest = args.iter();
        while let Some(&a) = rest.next() {
            let flag = |v: Option<&&str>| v.is_some_and(|v| *v != "0");
            match a {
                "-e" => self.simplify.0 = flag(rest.next()),
                "-f" => self.simplify.1 = flag(rest.next()),
                "-a" => {
                    rest.next();
                }
                _ => return Err(DrawError::Failed(format!("bsimplify: unknown option {a}"))),
            }
        }
        Ok(String::new())
    }
    /// DRAW command `bsmooth`: 1 uses in 1 converted scripts.
    ///
    /// Example: `bsmooth aSmoothCurve 1.0 …`
    pub fn bsmooth(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("bsmooth", args)
    }
    /// DRAW command `bsplinecurve`: 70 uses in 69 converted scripts.
    ///
    /// Example: `bsplinecurve w1 3 2 -1.0 4 1.0 4 0 0 0 1 2 10 0 1 4 6 0 1 10 0 0 1`
    pub fn bsplinecurve(&mut self, args: &[&str]) -> DrawResult {
        // `bsplinecurve name degree nknots k1 m1 ... x1 y1 z1 w1 ...`.
        need(args, 3, "bsplinecurve name degree nknots knots... poles...")?;
        let degree = self.number(args[1])? as usize;
        let nk = self.number(args[2])? as usize;
        let v = self.numbers(&args[3..])?;
        if v.len() < 2 * nk {
            return Err(DrawError::Failed("bsplinecurve: missing knots".into()));
        }
        let knots: Vec<f64> = (0..nk).map(|i| v[2 * i]).collect();
        let mults: Vec<usize> = (0..nk).map(|i| v[2 * i + 1] as usize).collect();
        let rest = &v[2 * nk..];
        if rest.len() % 4 != 0 {
            return Err(DrawError::Failed("bsplinecurve: poles come as x y z w".into()));
        }
        let poles: Vec<Point> = rest.chunks(4).map(|c| Point::new(c[0], c[1], c[2])).collect();
        let weights: Vec<f64> = rest.chunks(4).map(|c| c[3]).collect();
        let mut curve = BSplineCurve::new(poles, knots, mults, degree);
        if weights.iter().any(|&w| w != 1.0) {
            curve = curve.with_weights(weights);
        }
        self.set(args[0], DrawObject::Curve(Curve::BSpline(curve)));
        Ok(String::new())
    }
    /// DRAW command `bsplinesurf`: 7 uses in 7 converted scripts.
    ///
    /// Example: `bsplinesurf s 2 4 0.0 3 0.34 1 0.67 1 1.0 3 1 2 0.0 2 1.0 2 10 0 0 1 6 0 0 1 4 0 0 1 2 0 0 1 0 0 0 1`
    pub fn bsplinesurf(&mut self, args: &[&str]) -> DrawResult {
        // `bsplinesurf name udeg nuk uk um ... vdeg nvk vk vm ... poles (x y z w), u varying
        // fastest`.
        need(args, 3, "bsplinesurf name udeg nuk knots... vdeg nvk knots... poles...")?;
        let v = self.numbers(&args[1..])?;
        let knots = |at: usize| -> DrawResult<(usize, Vec<f64>, Vec<usize>, usize)> {
            let (degree, n) = (*v.get(at).unwrap_or(&0.0) as usize, *v.get(at + 1).unwrap_or(&0.0) as usize);
            let list = v.get(at + 2..at + 2 + 2 * n).ok_or_else(|| DrawError::Failed("bsplinesurf: missing knots".into()))?;
            Ok((degree, list.chunks(2).map(|c| c[0]).collect(), list.chunks(2).map(|c| c[1] as usize).collect(), at + 2 + 2 * n))
        };
        let (udeg, uk, um, at) = knots(0)?;
        let (vdeg, vk, vm, at) = knots(at)?;
        let nu = um.iter().sum::<usize>().saturating_sub(udeg + 1);
        let nv = vm.iter().sum::<usize>().saturating_sub(vdeg + 1);
        let rest = &v[at..];
        if rest.len() != 4 * nu * nv {
            return Err(DrawError::Failed("bsplinesurf: wrong number of poles".into()));
        }
        let pole = |i: usize, j: usize| &rest[(j * nu + i) * 4..][..4];
        let poles = (0..nu).map(|i| (0..nv).map(|j| { let c = pole(i, j); Point::new(c[0], c[1], c[2]) }).collect()).collect();
        let weights: Vec<Vec<f64>> = (0..nu).map(|i| (0..nv).map(|j| pole(i, j)[3]).collect()).collect();
        let mut surface = BSplineSurface::try_new(poles, KnotVector::new(uk, um, udeg), KnotVector::new(vk, vm, vdeg))?;
        if weights.iter().flatten().any(|&w| w != 1.0) {
            surface = surface.with_weights(weights);
        }
        self.set(args[0], DrawObject::Surface(Surface::BSpline(surface)));
        Ok(String::new())
    }
    /// DRAW command `bsplit`: 13 uses in 13 converted scripts.
    ///
    /// Example: `bsplit result`
    pub fn bsplit(&mut self, args: &[&str]) -> DrawResult {
        // `bsplit result [-t]`: the objects split by one another and by the tools: their
        // pieces.
        need(args, 1, "bsplit result")?;
        let objects = self.fuse_objects.len();
        let pieces: Vec<Shape> = self.fuse_cells()?.iter().filter(|c| c.inside[..objects].iter().any(|&i| i)).map(|c| c.shape).collect();
        let result = self.compound_of(&pieces)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `btolx`: 16 uses in 10 converted scripts.
    ///
    /// Example: `btolx tool`
    pub fn btolx(&mut self, args: &[&str]) -> DrawResult {
        // `btolx shape`: its tolerances corrected to cover its geometry's gaps.
        need(args, 1, "btolx shape")?;
        let shape = self.shape_var(args[0])?;
        modify::update_tolerances(&mut self.model, shape)?;
        Ok(String::new())
    }
    /// DRAW command `build3d`: 2 uses in 2 converted scripts.
    ///
    /// Example: `build3d s1`
    pub fn build3d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("build3d", args)
    }
    /// DRAW command `buildevol`: 103 uses in 103 converted scripts.
    ///
    /// Example: `buildevol`
    pub fn buildevol(&mut self, args: &[&str]) -> DrawResult {
        // `buildevol`: builds the evolving fillet described by mkevol and updatevol.
        let Some((name, shape, edges)) = self.evolving.take() else {
            return Err(DrawError::Failed("buildevol: no mkevol before".into()));
        };
        let _ = args;
        let result = blend::fillet_evolving(&mut self.model, shape, &edges)?;
        self.set_shape_var(&name, result)
    }
    /// DRAW command `buildsweep`: 108 uses in 108 converted scripts.
    ///
    /// Example: `buildsweep result -C -S`
    pub fn buildsweep(&mut self, args: &[&str]) -> DrawResult {
        // `buildsweep result [-M|-C|-R] [-S]`: the profile swept; a solid with -S (closed
        // by its ends), the swept shell otherwise. At corners of the spine the frame is
        // kept (-M, OCCT's default), or turned with right (-C) or round (-R) corners. As in
        // DRAW, only the argument after the result names a transition, and -S only follows.
        need(args, 1, "buildsweep result [-M|-C|-R] [-S]")?;
        let mut next = 1;
        let corner = match args.get(1).copied() {
            Some("-C") => kentron::sweep::Corner::Mitred,
            Some("-R") => kentron::sweep::Corner::Rounded,
            Some("-M") => kentron::sweep::Corner::Unturned,
            _ => {
                next = 0;
                kentron::sweep::Corner::Unturned
            }
        };
        let solid = args.get(1 + next) == Some(&"-S");
        let Some(spine) = self.sweep_spine else {
            return self.not_implemented("buildsweep", args);
        };
        if self.sweep_profiles.is_empty() {
            return Err(DrawError::Failed("buildsweep: no sections".into()));
        }
        let corners = {
            let mut tangents = Vec::new();
            for e in self.model.explore(spine, NodeKind::Edge).collect::<Vec<_>>() {
                if let Some((c, a, b)) = self.model.edge_curve(e) {
                    let (a, b) = if e.orientation == kentron::topology::Orientation::Reversed { (b, a) } else { (a, b) };
                    let sign = if b < a { -1.0 } else { 1.0 };
                    tangents.push((c.d1(a).1 * sign, c.d1(b).1 * sign));
                }
            }
            tangents.windows(2).any(|w| {
                let (u, v) = (w[0].1, w[1].0);
                u.dot(v) < (1.0 - 1e-6) * u.length() * v.length()
            })
        };
        if corners && corner != kentron::sweep::Corner::Unturned {
            let [(profile, contact, correction)] = self.sweep_profiles[..] else {
                return self.not_implemented("buildsweep", args);
            };
            // (Frenet's frame turns as the corrected one does where the spine is straight
            // between its corners.)
            let straight = self.model.explore(spine, NodeKind::Edge).all(|e| matches!(self.model.edge_curve(e), Some((kentron::geometry::curve::Curve::Line(_), _, _))));
            if !self.sweep_mode.is_empty() && self.sweep_mode != ["-CF"] && !(self.sweep_mode == ["-FR"] && straight) {
                return self.not_implemented("buildsweep", args);
            }
            let profile = if solid && profile.kind() == NodeKind::Wire { make::planar_face(&mut self.model, profile)? } else { profile };
            let options = kentron::sweep::PipeOptions { contact, correction, ..Default::default() };
            let r = kentron::sweep::pipe_cornered(&mut self.model, profile, spine, &options, corner)?;
            return self.set_shape_var(args[0], r);
        }
        if self.sweep_mode == ["-FR"] {
            // Frenet's frame only where it is the rotation-minimizing one: on a plane.
            let wire = if spine.kind() == NodeKind::Wire { spine } else { make::wire(&mut self.model, &[spine])? };
            if make::plane_of_wire(&self.model, wire).is_err() {
                return self.not_implemented("buildsweep", args);
            }
        }
        let r = if let [(profile, contact, correction)] = self.sweep_profiles[..] {
            let profile = if solid && profile.kind() == NodeKind::Wire { make::planar_face(&mut self.model, profile)? } else { profile };
            let options = kentron::sweep::PipeOptions { contact, correction, ..Default::default() };
            kentron::sweep::pipe_with(&mut self.model, profile, spine, &options)?.shape
        } else {
            let sections: Vec<(Shape, kentron::sweep::PipeOptions)> =
                self.sweep_profiles.iter().map(|&(s, contact, correction)| (s, kentron::sweep::PipeOptions { contact, correction, ..Default::default() })).collect();
            kentron::sweep::pipe_through(&mut self.model, &sections, spine, solid)?
        };
        self.set_shape_var(args[0], r)
    }
    /// DRAW command `isbbinterf`.
    ///
    /// Example: `isbbinterf b1 b2 -o`
    pub fn isbbinterf(&mut self, args: &[&str]) -> DrawResult {
        // `isbbinterf shape1 shape2 [-o]`: whether their bounding boxes (oriented ones with
        // `-o`) overlap.
        need(args, 2, "isbbinterf shape1 shape2 [-o]")?;
        let (a, b) = (self.shape_var(args[0])?, self.shape_var(args[1])?);
        if args.get(2) == Some(&"-o") {
            let (x, y) = (properties::oriented_bounding_box(&self.model, a), properties::oriented_bounding_box(&self.model, b));
            return Ok(if x.intersects(&y) { "The shapes are interfered by OBB.\n" } else { "The shapes are NOT interfered by OBB.\n" }.into());
        }
        let (x, y) = (properties::bounding_box(&self.model, a), properties::bounding_box(&self.model, b));
        let apart = x.max.x < y.min.x || y.max.x < x.min.x || x.max.y < y.min.y || y.max.y < x.min.y || x.max.z < y.min.z || y.max.z < x.min.z;
        Ok(if apart { "The shapes are NOT interfered by AABB.\n" } else { "The shapes are interfered by AABB.\n" }.into())
    }
    /// DRAW command `cmovetangent`.
    ///
    /// Example: `cmovetangent b1 2.0 x y dx dy 1`
    pub fn cmovetangent(&mut self, args: &[&str]) -> DrawResult {
        // `cmovetangent name u x y [z] tx ty [tz] [constraint]`: the B-spline moved to pass
        // through the point at `u` with the tangent there, keeping `constraint + 1` poles at
        // each end (one by default).
        need(args, 6, "cmovetangent name u x y [z] tx ty [tz] [constraint]")?;
        let v = self.numbers(&args[1..])?;
        let tolerance = 1e-5;
        let constraint = |c: Option<&f64>, degree: usize| c.map_or(0, |&c| (c as i32).max(-1).min(degree as i32 - 1));
        let moved = if args.len() <= 7 {
            let DrawObject::Curve2(Curve2::BSpline(mut b)) = self.get(args[0]).cloned().ok_or_else(|| DrawError::Failed(format!("{} is not a DRAW variable", args[0])))? else {
                return Ok(String::new());
            };
            let c = constraint(v.get(5), b.degree());
            let ok = b.move_point_and_tangent(v[0], Point2::new(v[1], v[2]), Vector2::new(v[3], v[4]), tolerance, c, c);
            self.set(args[0], DrawObject::Curve2(Curve2::BSpline(b)));
            ok
        } else {
            need(args, 8, "cmovetangent name u x y z tx ty tz [constraint]")?;
            let DrawObject::Curve(Curve::BSpline(mut b)) = self.get(args[0]).cloned().ok_or_else(|| DrawError::Failed(format!("{} is not a DRAW variable", args[0])))? else {
                return Ok(String::new());
            };
            let c = constraint(v.get(7), b.degree());
            let ok = b.move_point_and_tangent(v[0], Point::new(v[1], v[2], v[3]), Vector::new(v[4], v[5], v[6]), tolerance, c, c);
            self.set(args[0], DrawObject::Curve(Curve::BSpline(b)));
            ok
        };
        Ok(if moved.is_some() { String::new() } else { "Not enough degree of freedom increase degree please\n".into() })
    }
    /// DRAW command `setaxis`.
    ///
    /// Example: `setaxis 0 0 0 0 0 1 1 0 0`
    pub fn setaxis(&mut self, args: &[&str]) -> DrawResult {
        // `setaxis x y z nx ny nz xx xy xz`: the axis (and x direction) of the helices made.
        need(args, 9, "setaxis x y z nx ny nz xx xy xz")?;
        let v = self.numbers(&args[..9])?;
        self.helix_frame = Some(Self::curve_frame(&v)?);
        Ok(String::new())
    }
    /// DRAW command `helix`.
    ///
    /// Example: `helix name np D1 H1 [Hi...] P1 [Pi...] PF1 [PFi...]`
    pub fn helix(&mut self, args: &[&str]) -> DrawResult {
        // `helix name np D1 H1 [Hi...] P1 [Pi...] PF1 [PFi...]`.
        self.helix_command("helix", args)
    }
    /// DRAW command `helix2`.
    ///
    /// Example: `helix2 name np D1 P1 [Pi...] N1 [Ni...]`
    pub fn helix2(&mut self, args: &[&str]) -> DrawResult {
        // `helix2 name np D1 P1 [Pi...] N1 [Ni...]`.
        self.helix_command("helix2", args)
    }
    /// DRAW command `comphelix`.
    ///
    /// Example: `comphelix name np D1 [Di...] H1 [Hi...] P1 [Pi...] PF1 [PFi...]`
    pub fn comphelix(&mut self, args: &[&str]) -> DrawResult {
        // `comphelix name np D1 [Di...] H1 [Hi...] P1 [Pi...] PF1 [PFi...]`.
        self.helix_command("comphelix", args)
    }
    /// DRAW command `comphelix2`.
    ///
    /// Example: `comphelix2 name np D1 [Di...] P1 [Pi...] N1 [Ni...]`
    pub fn comphelix2(&mut self, args: &[&str]) -> DrawResult {
        // `comphelix2 name np D1 [Di...] P1 [Pi...] N1 [Ni...]`.
        self.helix_command("comphelix2", args)
    }
    /// DRAW command `spiral`.
    ///
    /// Example: `spiral name np D1 D2 H1 [Hi...] P1 [Pi...] PF1 [PFi...]`
    pub fn spiral(&mut self, args: &[&str]) -> DrawResult {
        // `spiral name np D1 D2 H1 [Hi...] P1 [Pi...] PF1 [PFi...]`.
        self.helix_command("spiral", args)
    }
    /// DRAW command `spiral2`.
    ///
    /// Example: `spiral2 name np D1 D2 P1 [Pi...] N1 [Ni...]`
    pub fn spiral2(&mut self, args: &[&str]) -> DrawResult {
        // `spiral2 name np D1 D2 P1 [Pi...] N1 [Ni...]`.
        self.helix_command("spiral2", args)
    }
    /// DRAW command `chamf_throat`.
    ///
    /// Example: `chamf_throat result b b_8 2.`
    pub fn chamf_throat(&mut self, args: &[&str]) -> DrawResult {
        // `chamf_throat result shape edge throat [edge throat ...]`: symmetric chamfers whose
        // faces lie at `throat` from their edges, between faces meeting at the angle at the
        // edge's middle (exactly so between planes).
        use kentron::blend::{ChamferSize, ChamferedEdge};
        use kentron::topology::{EdgeId, FaceId};
        need(args, 4, "chamf_throat result shape edge throat [edge throat ...]")?;
        let shape = self.shape_var(args[1])?;
        let faces: Vec<Shape> = self.model.explore(shape, NodeKind::Face).collect();
        let mut edges = Vec::new();
        for pair in args[2..].chunks(2) {
            let [edge, throat] = pair else { return Err(DrawError::Failed("chamf_throat: an edge without a throat".into())) };
            let (edge, throat) = (self.shape_var(edge)?, self.number(throat)?);
            let edge_id = EdgeId::try_from(edge.node).map_err(|_| DrawError::Failed("chamf_throat: not an edge".into()))?;
            // The outward normals of the faces along the edge, at its middle.
            let mut normals = Vec::new();
            for &face in &faces {
                let Ok(face_id) = FaceId::try_from(face.node) else { continue };
                let Some(coedge) = self.model.find_coedge(edge_id, face_id) else { continue };
                let (Some((pcurve, first, last)), Some(surface)) = (self.model.effective_pcurve(coedge), self.model.surface(face)) else { continue };
                let uv = pcurve.point(0.5 * (first + last));
                let Some(n) = kentron::geometry::props::surface_normal(surface, uv.x, uv.y, 1e-12) else { continue };
                let n = Vector::from(n.direction());
                normals.push(if face.orientation == Orientation::Reversed { -n } else { n });
            }
            let [a, b] = normals[..] else {
                return Err(DrawError::Failed("chamf_throat: the edge must bound two faces".into()));
            };
            // The faces' angle inside the solid is π less the normals'.
            let inside = std::f64::consts::PI - a.angle(b);
            let distance = throat / (0.5 * inside).cos();
            edges.push(ChamferedEdge { edge, size: ChamferSize::Distances(distance, distance), face: None });
        }
        let result = kentron::blend::chamfer_edges(&mut self.model, shape, &edges)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `chamf`: 17 uses in 8 converted scripts.
    ///
    /// Example: `chamf result s s_10_3 50`
    pub fn chamf(&mut self, args: &[&str]) -> DrawResult {
        // `chamf result shape` then edges, each `edge dist`, or `edge face` and `S dist`,
        // `dist1 dist2` or `A dist angle` (degrees), the first distance along the face.
        use kentron::blend::{ChamferSize, ChamferedEdge};
        if args.len() < 4 {
            return self.not_implemented("chamf", args);
        }
        let shape = self.shape_var(args[1])?;
        let mut edges = Vec::new();
        let mut i = 2;
        while i < args.len() {
            let edge = self.shape_var(args[i])?;
            i += 1;
            let rest = &args[i..];
            let (size, face, used) = if let Some(Ok(d)) = rest.first().map(|a| self.number(a)) {
                (ChamferSize::Distances(d, d), None, 1)
            } else {
                let face = Some(self.shape_var(rest.first().copied().unwrap_or(""))?);
                match rest.get(1).copied() {
                    Some("S") => {
                        let d = self.number(rest.get(2).copied().unwrap_or(""))?;
                        (ChamferSize::Distances(d, d), face, 3)
                    }
                    Some("A") => {
                        let d = self.number(rest.get(2).copied().unwrap_or(""))?;
                        let angle = self.number(rest.get(3).copied().unwrap_or(""))?;
                        (ChamferSize::DistanceAngle(d, angle.to_radians()), face, 4)
                    }
                    _ => {
                        let d1 = self.number(rest.get(1).copied().unwrap_or(""))?;
                        let d2 = self.number(rest.get(2).copied().unwrap_or(""))?;
                        (ChamferSize::Distances(d1, d2), face, 3)
                    }
                }
            };
            edges.push(ChamferedEdge { edge, size, face });
            i += used;
        }
        let result = kentron::blend::chamfer_edges(&mut self.model, shape, &edges)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `chamfer2d`: 1 uses in 1 converted scripts.
    ///
    /// Example: `chamfer2d result e1 e2 10 20`
    pub fn chamfer2d(&mut self, args: &[&str]) -> DrawResult {
        // `chamfer2d result edge1 edge2 d1 d2`: the chamfer segment (the edges cut back to it), or
        // `chamfer2d result wire d1 d2`: the wire of two edges with its corner bevelled.
        self.corner_2d("chamfer2d", args, 2)
    }
    /// DRAW command `checkcurveonsurf`: 1 uses in 1 converted scripts.
    ///
    /// Example: `checkcurveonsurf f`
    pub fn checkcurveonsurf(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("checkcurveonsurf", args)
    }
    /// DRAW command `checkdump`: 46 uses in 46 converted scripts.
    ///
    /// Example: `checkdump result -name Center Axis XAxis YAxis Radii -ref {0 0 0} {0 0 1} {1 0 -0} {-0 1 0} {20 10} `
    pub fn checkdump(&mut self, args: &[&str]) -> DrawResult {
        // `checkdump name -name {fields} -ref {{values} ...} [-eps e]`: compares fields of a
        // line's or conic's description (Origin, Center, Axis, XAxis, YAxis, Radii, Focal)
        // as CheckCommands.tcl does: relatively within `e`, else exactly.
        need(args, 1, "checkdump name -name fields -ref values [-eps e]")?;
        let option = |name: &str| args.iter().position(|a| *a == name).and_then(|i| args.get(i + 1)).copied();
        let names: Vec<&str> = option("-name").unwrap_or("").split_whitespace().collect();
        let refs = option("-ref").unwrap_or("");
        let groups: Vec<Vec<f64>> = if names.len() > 1 {
            refs.split('}').map(|g| g.trim().trim_start_matches('{').split_whitespace().filter_map(|x| x.parse().ok()).collect::<Vec<f64>>()).filter(|g| !g.is_empty()).collect()
        } else {
            vec![refs.split_whitespace().filter_map(|x| x.trim_matches(|c| c == '{' || c == '}').parse().ok()).collect()]
        };
        let epsilon = option("-eps").map(|e| self.number(e)).transpose()?;
        let fields = self.dump_fields(args[0])?;
        let mut errors = Vec::new();
        for (name, expected) in names.iter().zip(&groups) {
            let Some(actual) = fields.iter().find(|f| f.0 == *name).map(|f| &f.1) else {
                errors.push(format!("{name} is not listed in the dump"));
                continue;
            };
            for (&item, &res) in expected.iter().zip(actual) {
                let wrong = match epsilon {
                    None => item != res,
                    Some(eps) => (res.abs() > 1e-7 || item.abs() > 1e-7) && if item != 0.0 { ((item - res) / item).abs() > eps } else { res != 0.0 },
                };
                if wrong {
                    errors.push(format!("{name}: {res} instead of {item}"));
                }
            }
            if actual.len() < expected.len() {
                errors.push(format!("{name}: {actual:?} has too few values"));
            }
        }
        match errors.is_empty() {
            true => Ok(String::new()),
            false => Err(DrawError::Failed(format!("Error: checkdump {}: {}", args[0], errors.join("; ")))),
        }
    }
    /// DRAW command `checkfaults`: 700 uses in 700 converted scripts.
    ///
    /// Example: `checkfaults result a 0`
    pub fn checkfaults(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("checkfaults", args)
    }
    /// DRAW command `checkfreebounds`: 739 uses in 738 converted scripts.
    ///
    /// Example: `checkfreebounds result 2`
    pub fn checkfreebounds(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("checkfreebounds", args)
    }
    /// DRAW command `checklength`: 181 uses in 174 converted scripts.
    ///
    /// Example: `checklength o -l 2.3717833300483151`
    pub fn checklength(&mut self, args: &[&str]) -> DrawResult {
        // `checklength name -l value [-eps e]`: the length of a curve (over its domain) or of
        // a shape's edges, within the relative tolerance 1e-4 (CheckCommands.tcl's default).
        need(args, 3, "checklength name -l value")?;
        if args[1] != "-l" {
            return self.not_implemented("checklength", args);
        }
        let expected = self.number(args[2])?;
        let relative = match args.iter().position(|a| *a == "-eps") {
            Some(i) => self.number(args.get(i + 1).copied().unwrap_or("1e-4"))?,
            None => 1e-4,
        };
        let actual = match self.get(args[0]).cloned() {
            Some(DrawObject::Shape(s)) => properties::length(&self.model, s).mass,
            Some(DrawObject::Curve2(c)) => kentron::geometry::measure::length2(&c),
            _ => {
                let curve = self.curve_var(args[0])?;
                let edge = make::edge_from_curve(&mut self.model, curve)?;
                properties::length(&self.model, edge).mass
            }
        };
        match (actual - expected).abs() <= relative * expected.abs().max(1e-300) {
            true => Ok(String::new()),
            false => Err(DrawError::Failed(format!("Error: The length of {} is {actual}, expected {expected}", args[0]))),
        }
    }
    /// DRAW command `checkmaxtol`: 766 uses in 760 converted scripts.
    ///
    /// Example: `checkmaxtol result -ref 0.001`
    pub fn checkmaxtol(&mut self, args: &[&str]) -> DrawResult {
        // `checkmaxtol shape [-ref value] [-min_tol t] [-source shapes] [-multi_tol k]` (OCCT's
        // CheckCommands.tcl): the largest tolerance of the shape's sub-shapes, an error if above
        // the checking tolerance or off the reference.
        need(args, 1, "checkmaxtol shape [options]")?;
        let max_of = |d: &Self, s: Shape| -> f64 {
            [NodeKind::Vertex, NodeKind::Edge, NodeKind::Face]
                .into_iter()
                .flat_map(|k| d.model.unique_subshapes(s, k))
                .map(|x| d.model.tolerance(x))
                .fold(0.0, f64::max)
        };
        let shape = self.shape_var(args[0])?;
        let max = max_of(self, shape);
        let (mut reference, mut min_tol, mut multiplier, mut sources) = (None, 0.0, 0.0, Vec::new());
        let mut k = 1;
        while k < args.len() {
            let value = args.get(k + 1).copied().ok_or_else(|| DrawError::Failed(format!("checkmaxtol: {} needs a value", args[k])))?;
            match args[k] {
                "-ref" => reference = Some(self.number(value)?),
                "-min_tol" => min_tol = self.number(value)?,
                "-multi_tol" => multiplier = self.number(value)?,
                "-source" => sources = value.split_whitespace().map(|n| self.shape_var(n)).collect::<DrawResult<Vec<_>>>()?,
                other => return Err(DrawError::Failed(format!("Error: unsupported option \"{other}\""))),
            }
            k += 2;
        }
        let mut source_max = 0.0f64;
        for s in sources {
            source_max = max_of(self, s);
            min_tol = f64::max(min_tol, source_max);
        }
        if multiplier != 0.0 {
            min_tol = multiplier * source_max;
        }
        if min_tol != 0.0 && max > min_tol {
            return Err(DrawError::Failed(format!("Error: tolerance of \"{}\" ({}) is greater than checking tolerance ({})", args[0], format_real(max), format_real(min_tol))));
        }
        if let Some(r) = reference {
            self.checkreal(&["Max tolerance", &format_real(max), &format_real(r), "0.0001", "0.01"])?;
        }
        Ok(format_real(max))
    }
    /// DRAW command `checknbshapes`: 2930 uses in 2542 converted scripts.
    ///
    /// Example: `checknbshapes result -vertex 4 -edge 2 -t`
    pub fn checknbshapes(&mut self, args: &[&str]) -> DrawResult {
        // `checknbshapes name -vertex n -edge n ... -shape n`: the numbers of distinct
        // sub-shapes of each kind (the shape itself included), as OCCT's `nbshapes` reports.
        need(args, 1, "checknbshapes name -kind count ...")?;
        let shape = self.shape_var(args[0])?;
        let kinds = SHAPE_KINDS.map(|(name, kind)| (name, kind));
        let shared = !args.contains(&"-t");
        let count = |kind: NodeKind| count_subshapes(&self.model, shape, kind, shared);
        let mut errors = Vec::new();
        let mut i = 1;
        while i < args.len() {
            if args[i] == "-t" {
                i += 1;
                continue;
            }
            // (`-m message`: what the check is about, for its report.)
            if args[i] == "-m" {
                i += 2;
                continue;
            }
            // (`-ref report`: the counts as `nbshapes` reports them, a `KIND : n` line each;
            // lines without a number are passed over, as OCCT's pattern does.)
            if args[i] == "-ref" {
                let report = args.get(i + 1).ok_or_else(|| DrawError::Failed("-ref needs a report".into()))?;
                for line in report.lines() {
                    let Some((name, value)) = line.split_once(':') else { continue };
                    let Ok(expected) = value.trim().parse::<usize>() else { continue };
                    let name = name.trim();
                    let actual = match name {
                        "SHAPE" => kinds.iter().map(|&(_, k)| count(k)).sum(),
                        _ => match kinds.iter().find(|(o, _)| o.eq_ignore_ascii_case(name)) {
                            Some(&(_, kind)) => count(kind),
                            None => continue,
                        },
                    };
                    if actual != expected {
                        errors.push(format!("{}: {actual} instead of {expected}", name.to_lowercase()));
                    }
                }
                i += 2;
                continue;
            }
            let expected = args
                .get(i + 1)
                .ok_or_else(|| DrawError::Failed(format!("{} needs a count", args[i])))?
                .parse::<usize>()
                .map_err(|_| DrawError::Failed(format!("{} needs a count", args[i])))?;
            let actual = match args[i] {
                "-shape" => kinds.iter().map(|&(_, k)| count(k)).sum(),
                option => match kinds.iter().find(|(o, _)| option.strip_prefix('-') == Some(*o)) {
                    Some(&(_, kind)) => count(kind),
                    None => return self.not_implemented("checknbshapes", args),
                },
            };
            if actual != expected {
                errors.push(format!("{}: {actual} instead of {expected}", &args[i][1..]));
            }
            i += 2;
        }
        match errors.is_empty() {
            true => Ok(String::new()),
            false => Err(DrawError::Failed(format!("Error : wrong numbers of sub-shapes of {}: {}", args[0], errors.join(", ")))),
        }
    }
    /// DRAW command `checkprops`: 6565 uses in 6078 converted scripts.
    ///
    /// Example: `checkprops result -s 59527.9`
    pub fn checkprops(&mut self, args: &[&str]) -> DrawResult {
        // `checkprops shape [-l|-s|-v value] [-eps e] [-deps d]`: compares the length, area or
        // volume with the expected value within the relative tolerance `-deps` (1e-2), as
        // OCCT's CheckCommands.tcl does; `empty` expects exactly zero. `-eps` is a computation
        // precision and is ignored.
        need(args, 1, "checkprops shape [-l|-s|-v value] [-eps e] [-deps d]")?;
        let shape = self.shape_var(args[0])?;
        let mut checks = Vec::new();
        let mut relative = 1.0e-2;
        let mut i = 1;
        while i < args.len() {
            let value = || args.get(i + 1).ok_or_else(|| DrawError::Failed(format!("{} needs a value", args[i])));
            match args[i] {
                "-l" | "-s" | "-v" => match *value()? {
                    "empty" => checks.push((args[i], None)),
                    v => checks.push((args[i], Some(self.number(v)?))),
                },
                "-deps" => relative = self.number(value()?)?,
                "-eps" => {}
                _ => return self.not_implemented("checkprops", args),
            }
            i += 2;
        }
        let mut out = String::new();
        for (kind, expected) in checks {
            let (name, actual) = match kind {
                "-l" => ("length", properties::length(&self.model, shape).mass),
                "-s" => ("area", properties::area(&self.model, shape).mass),
                "-v" => ("volume", properties::volume(&self.model, shape).mass),
                _ => return self.not_implemented("checkprops", args),
            };
            match expected {
                None if actual != 0.0 => {
                    return Err(DrawError::Failed(format!("Error : The command is not valid. The {name} is {actual}")));
                }
                Some(_) if actual == 0.0 => {
                    return Err(DrawError::Failed(format!("Error : The command is not valid. The {name} is 0.")));
                }
                Some(expected)
                    if (expected != 0.0 && ((expected - actual) / expected).abs() > relative)
                        || (expected == 0.0 && actual != 0.0) =>
                {
                    return Err(DrawError::Failed(format!(
                        "Error : The {name} of result shape is {actual}, expected {expected}"
                    )));
                }
                _ => {}
            }
            out += &format!("{name}: {actual}\n");
        }
        Ok(out)
    }
    /// DRAW command `checkreal`: 59 uses in 20 converted scripts.
    ///
    /// Example: `checkreal Distance 1 … … … …`
    pub fn checkreal(&mut self, args: &[&str]) -> DrawResult {
        // `checkreal name value expected tol_abs tol_rel` (OCCT's CheckCommands.tcl): an error
        // unless the value is within either tolerance of the expected one.
        need(args, 5, "checkreal name value expected tol_abs tol_rel")?;
        let (value, expected) = (self.number(args[1])?, self.number(args[2])?);
        let (abs, rel) = (self.number(args[3])?, self.number(args[4])?);
        let diff = (value - expected).abs();
        if diff > abs && diff > rel * expected.abs() {
            return Err(DrawError::Failed(format!("Error: {} = {} is not equal to expected {}", args[0], format_real(value), format_real(expected))));
        }
        Ok(format!("Check of {} OK", args[0]))
    }
    /// DRAW command `checksection`: 1095 uses in 1085 converted scripts.
    ///
    /// Example: `checksection result`
    pub fn checksection(&mut self, args: &[&str]) -> DrawResult {
        // `checksection shape [-r count]`: the number of vertices on only one edge end ("alone
        // vertices"), checked against `count` where given (otherwise only reported, as OCCT
        // does).
        need(args, 1, "checksection shape [-r count]")?;
        let shape = self.shape_var(args[0])?;
        let expected = match args.iter().position(|a| *a == "-r") {
            Some(i) => Some(self.number(args.get(i + 1).copied().unwrap_or("0"))? as usize),
            None => None,
        };
        let mut ends: Vec<(Shape, usize)> = Vec::new();
        for edge in self.model.unique_subshapes(shape, NodeKind::Edge) {
            let (a, b) = self.model.edge_vertices(edge);
            for v in [a, b] {
                match ends.iter_mut().find(|(w, _)| w.is_same(&v)) {
                    Some(e) => e.1 += 1,
                    None => ends.push((v, 1)),
                }
            }
        }
        let alone = ends.iter().filter(|e| e.1 == 1).count();
        match expected {
            Some(expected) if alone != expected => Err(DrawError::Failed(format!("Error: checksection {}: {alone} alone vertices, expected {expected}", args[0]))),
            _ => Ok(format!(" nb alone Vertices : {alone}")),
        }
    }
    /// DRAW command `checkshape`: 5472 uses in 4850 converted scripts.
    ///
    /// Example: `checkshape result`
    pub fn checkshape(&mut self, args: &[&str]) -> DrawResult {
        // `checkshape name [options]`: fails unless the shape is valid.
        // `-top` (check the shape itself, not its sub-shapes separately) changes nothing here:
        // the whole shape is checked.
        let args: Vec<&str> = args.iter().copied().filter(|a| *a != "-top").collect();
        need(&args, 1, "checkshape name")?;
        // (The check OCCT's test scripts end with, of `result` where there is one.)
        if args[0] == "result" && self.get("result").is_none() {
            return Ok(String::new());
        }
        let shape = self.shape_var(args[0])?;
        let report = validity::check(&self.model, shape);
        match report.is_valid() {
            true => Ok("This shape seems to be valid".into()),
            false => {
                let problems: Vec<String> = report.errors().map(|p| p.description.clone()).collect();
                Err(DrawError::Failed(format!("Faulty shape {}: {}", args[0], problems.join("; "))))
            }
        }
    }
    /// DRAW command `checktrinfo`: 186 uses in 173 converted scripts.
    ///
    /// Example: `checktrinfo result -tri 294 -nod 300`
    pub fn checktrinfo(&mut self, args: &[&str]) -> DrawResult {
        // `checktrinfo shape [-face n] [-tri [n]] [-nod [n]] [-defl [d]] [-tol_abs_* t]
        // [-tol_rel_* t]`: the number of faces, the numbers of triangles and nodes and the
        // largest deflection of the faces' meshes;
        // an option without a value checks that the number is not zero.
        need(args, 1, "checktrinfo shape [-tri [n]] [-nod [n]] [-defl [d]]")?;
        let shape = self.shape_var(args[0])?;
        let (mut triangles, mut nodes, mut deflection) = (0usize, 0usize, 0.0f64);
        let faces = self.model.unique_subshapes(shape, NodeKind::Face);
        for &face in &faces {
            if let Some(mesh) = self.model.face_mesh(face) {
                triangles += mesh.triangles().len();
                nodes += mesh.nodes().len();
                deflection = deflection.max(mesh.deflection());
            }
        }
        let mut tolerances = std::collections::HashMap::new();
        let mut checks: Vec<(&str, Option<f64>)> = Vec::new();
        let mut i = 1;
        while i < args.len() {
            let value = args.get(i + 1).filter(|v| !v.starts_with('-') || v.parse::<f64>().is_ok()).map(|v| self.number(v)).transpose()?;
            match args[i] {
                "-face" | "-tri" | "-nod" | "-defl" => checks.push((&args[i][1..], value)),
                t if t.starts_with("-tol_") => {
                    tolerances.insert(&t[5..], value.unwrap_or(0.0));
                }
                "-empty" => checks.push(("empty", None)),
                _ => return self.not_implemented("checktrinfo", args),
            }
            i += if value.is_some() { 2 } else { 1 };
        }
        let mut errors = Vec::new();
        for (what, expected) in checks {
            let actual = match what {
                "face" => faces.len() as f64,
                "tri" => triangles as f64,
                "nod" => nodes as f64,
                "defl" => deflection,
                _ => {
                    if triangles != 0 {
                        errors.push(format!("{triangles} triangles where none were expected"));
                    }
                    continue;
                }
            };
            let ok = match expected {
                None => actual != 0.0,
                Some(e) => {
                    let absolute = tolerances.get(format!("abs_{what}").as_str()).copied().unwrap_or(0.0);
                    let relative = tolerances.get(format!("rel_{what}").as_str()).copied().unwrap_or(0.0);
                    (actual - e).abs() <= absolute.max(relative * e.abs())
                }
            };
            if !ok {
                errors.push(format!("{what} is {actual}, expected {}", expected.map_or("non-zero".into(), |e| e.to_string())));
            }
        }
        match errors.is_empty() {
            true => Ok(String::new()),
            false => Err(DrawError::Failed(format!("Error: checktrinfo {}: {}", args[0], errors.join("; ")))),
        }
    }
    /// DRAW command `circ2d3Tan`: 3 uses in 3 converted scripts.
    ///
    /// Example: `circ2d3Tan tan1 c1 c2 l1`
    pub fn circ2d3Tan(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("circ2d3Tan", args)
    }
    /// DRAW command `circle`: 224 uses in 155 converted scripts.
    ///
    /// Example: `circle w -20 0 0 20`
    pub fn circle(&mut self, args: &[&str]) -> DrawResult {
        // `circle name x y [dx dy] r` (2D) or `circle name x y z [nx ny nz [ux uy uz]] r`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(1));
        let object = match frame_args.len() {
            2 | 4 => DrawObject::Curve2(Curve2::Circle(Circle2::new(Draw::curve_frame2(frame_args)?, r[0]))),
            _ => DrawObject::Curve(Curve::Circle(Circle::new(Draw::curve_frame(frame_args)?, r[0]))),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `cirtang`: 2 uses in 2 converted scripts.
    ///
    /// Example: `cirtang res -c c1 -c c2 -r .1`
    pub fn cirtang(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("cirtang", args)
    }
    /// DRAW command `clknots`: 1 uses in 1 converted scripts.
    ///
    /// Example: `clknots res`
    pub fn clknots(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("clknots", args)
    }
    /// DRAW command `clpoles`: 1 uses in 1 converted scripts.
    ///
    /// Example: `clpoles res`
    pub fn clpoles(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("clpoles", args)
    }
    /// DRAW command `compound`: 92 uses in 72 converted scripts.
    ///
    /// Example: `compound m_2 m_3 m_4 m_5 c1`
    pub fn compound(&mut self, args: &[&str]) -> DrawResult {
        // `compound [name ...] result`.
        need(args, 1, "compound [name ...] result")?;
        let n = args.len();
        // (Names holding no shape are passed over, as DRAW does.)
        let parts: Vec<Shape> = args[..n - 1].iter().filter(|a| self.get(a).is_some()).map(|a| self.shape_var(a)).collect::<DrawResult<Vec<Shape>>>()?;
        let compound = self.model.new_empty(NodeKind::Compound);
        for part in parts {
            self.model.add(compound, part)?;
        }
        self.set_shape_var(args[n - 1], compound)
    }
    /// DRAW command `computetolerance`: 4 uses in 3 converted scripts.
    ///
    /// Example: `computetolerance res1`
    pub fn computetolerance(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("computetolerance", args)
    }
    /// DRAW command `concatC0wire`: 7 uses in 4 converted scripts.
    ///
    /// Example: `concatC0wire aa a`
    pub fn concatC0wire(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("concatC0wire", args)
    }
    /// DRAW command `concatwire`: 1 uses in 1 converted scripts.
    ///
    /// Example: `concatwire r s`
    pub fn concatwire(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("concatwire", args)
    }
    /// DRAW command `cone`: 47 uses in 36 converted scripts.
    ///
    /// Example: `cone con_f2 0 112.0915332 -696.20434566999995 0 -0.32018607430384038 0.94735467372146087 79.81851702`
    pub fn cone(&mut self, args: &[&str]) -> DrawResult {
        // `cone name [x y z [nx ny nz [ux uy uz]]] semi-angle radius` (the angle in degrees).
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(2));
        let frame = if frame_args.is_empty() { Frame::WORLD } else { Draw::curve_frame(frame_args)? };
        self.set(args[0], DrawObject::Surface(Surface::Cone(kentron::geometry::surface::Cone::new(frame, r[0].to_radians(), r[1]))));
        Ok(String::new())
    }
    /// DRAW command `connectedges`: 1 uses in 1 converted scripts.
    ///
    /// Example: `connectedges result a`
    pub fn connectedges(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("connectedges", args)
    }
    /// DRAW command `convert`: 37 uses in 34 converted scripts.
    ///
    /// Example: `convert c c`
    pub fn convert(&mut self, args: &[&str]) -> DrawResult {
        // `convert name curve|surface`: the B-spline form of a curve, 2D curve or surface.
        need(args, 2, "convert name curve|surface")?;
        let failed = || DrawError::Failed(format!("convert: {} has no B-spline form", args[1]));
        let object = match self.get(args[1]).cloned() {
            Some(DrawObject::Curve(c)) => DrawObject::Curve(Curve::BSpline(gconvert::curve_to_bspline(&c).ok_or_else(failed)?)),
            Some(DrawObject::Curve2(c)) => DrawObject::Curve2(Curve2::BSpline(gconvert::curve2_to_bspline(&c).ok_or_else(failed)?)),
            Some(DrawObject::Surface(s)) => DrawObject::Surface(Surface::BSpline(gconvert::surface_to_bspline(&s).ok_or_else(failed)?)),
            _ => return self.not_implemented("convert", args),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `coordload`: 2 uses in 2 converted scripts.
    ///
    /// Example: `coordload res …`
    pub fn coordload(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("coordload", args)
    }
    /// DRAW command `copy`: 515 uses in 162 converted scripts.
    ///
    /// Example: `copy res_1 result`
    pub fn copy(&mut self, args: &[&str]) -> DrawResult {
        // `copy name newname [name newname ...]`: the variable's value under a new name (a
        // shape is shared, not duplicated: that is `tcopy`).
        if args.is_empty() || args.len() % 2 != 0 {
            return Err(DrawError::Failed("usage: copy name newname [name newname ...]".into()));
        }
        for pair in args.chunks(2) {
            let value = self.get(pair[0]).cloned().ok_or_else(|| DrawError::Failed(format!("{} is not a DRAW variable", pair[0])))?;
            self.set(pair[1], value);
        }
        Ok(String::new())
    }
    /// DRAW command `correctnormals`: 1 uses in 1 converted scripts.
    ///
    /// Example: `correctnormals a`
    pub fn correctnormals(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("correctnormals", args)
    }
    /// DRAW command `cprj`: 1 uses in 1 converted scripts.
    ///
    /// Example: `cprj res a f 100 100 20`
    pub fn cprj(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("cprj", args)
    }
    /// DRAW command `cvalue`: 5 uses in 5 converted scripts.
    ///
    /// Example: `cvalue c -0.58758239134636203 x y z`
    pub fn cvalue(&mut self, args: &[&str]) -> DrawResult {
        // `cvalue curve u x y z [dx dy dz [d2x d2y d2z]]`: the point (and derivatives) at u,
        // into numeric variables.
        need(args, 5, "cvalue curve u x y z [dx dy dz [d2x d2y d2z]]")?;
        let curve = self.curve_var(args[0])?;
        let u = self.number(args[1])?;
        let (p, d1, d2) = curve.d2(u);
        let values = [p.x, p.y, p.z, d1.x, d1.y, d1.z, d2.x, d2.y, d2.z];
        for (name, value) in args[2..].iter().zip(values) {
            self.set(name, DrawObject::Number(value));
        }
        Ok(String::new())
    }
    /// DRAW command `cylinder`: 79 uses in 53 converted scripts.
    ///
    /// Example: `cylinder cyl_f4 0 0 -4.829470157119431e-014 0 1 1.1102230246251565e-016 1750`
    pub fn cylinder(&mut self, args: &[&str]) -> DrawResult {
        // `cylinder name [x y z [nx ny nz [ux uy uz]]] radius`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(1));
        let frame = if frame_args.is_empty() { Frame::WORLD } else { Draw::curve_frame(frame_args)? };
        self.set(args[0], DrawObject::Surface(Surface::Cylinder(kentron::geometry::surface::Cylinder::new(frame, r[0]))));
        Ok(String::new())
    }
    /// DRAW command `dall`: 5 uses in 5 converted scripts.
    ///
    /// Example: `dall`
    pub fn dall(&mut self, _args: &[&str]) -> DrawResult {
        // Forgets every variable (the shapes stay in the model, unreferenced).
        self.vars.clear();
        Ok(String::new())
    }
    /// DRAW command `depouille`: 131 uses in 131 converted scripts.
    ///
    /// Example: `depouille result p 0 0 1 p_1 5 0 0 0 0 0 1 p_7 1 0 0 0 0 0 1 p_8 5 0 0 0 0 0 1 p_2 1 0 0 0 0 0 1`
    pub fn depouille(&mut self, args: &[&str]) -> DrawResult {
        // `depouille result shape dx dy dz face angle px py pz nx ny nz [face angle ...]`: the
        // pull direction, then per face its angle in degrees and its neutral plane.
        if args.len() < 13 || (args.len() - 5) % 8 != 0 {
            return self.not_implemented("depouille", args);
        }
        let shape = self.shape_var(args[1])?;
        let d = self.numbers(&args[2..5])?;
        let direction = Direction::try_from_vector(Vector::new(d[0], d[1], d[2]))
            .ok_or_else(|| DrawError::Failed("depouille: null direction".into()))?;
        let mut faces = Vec::new();
        for group in args[5..].chunks(8) {
            let face = self.shape_var(group[0])?;
            let v = self.numbers(&group[1..8])?;
            let normal = Direction::try_from_vector(Vector::new(v[4], v[5], v[6]))
                .ok_or_else(|| DrawError::Failed("depouille: null plane normal".into()))?;
            let neutral = kentron::geometry::surface::Plane::new(Point::new(v[1], v[2], v[3]), normal);
            faces.push(kentron::draft::DraftedFace { face, angle: v[0].to_radians(), neutral });
        }
        let result = kentron::draft::draft(&mut self.model, shape, direction, &faces)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `disp`: 1 uses in 1 converted scripts.
    ///
    /// Example: `disp l p`
    pub fn disp(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("disp", args)
    }
    /// DRAW command `distmini`: 21 uses in 16 converted scripts.
    ///
    /// Example: `distmini dv1 vertex1 ellipse`
    pub fn distmini(&mut self, args: &[&str]) -> DrawResult {
        // `distmini name shape1 shape2`: the distance as `name_val`, and as `name` the segment
        // between the closest points (a vertex if they meet).
        need(args, 3, "distmini name shape1 shape2")?;
        // (A name of no shape, as OCCT's scripts going through all DRAW variables give it:
        // reported, no error.)
        let (Ok(a), Ok(b)) = (self.shape_var(args[1]), self.shape_var(args[2])) else {
            return Ok("problem: no distance is found\n".into());
        };
        let d = kentron::distance::distance(&mut self.model, a, b)?;
        self.set(&format!("{}_val", args[0]), DrawObject::Number(d.value));
        let segment = if d.value > tolerance::LINEAR {
            make::edge_between(&mut self.model, d.on_first, d.on_second)?
        } else {
            make::vertex(&mut self.model, d.on_first)
        };
        self.set_shape_var(args[0], segment)?;
        Ok(format!("the distance value is : {}\nthe number of solutions is :1", format_real(d.value)))
    }
    /// DRAW command `draft`: 1 uses in 1 converted scripts.
    ///
    /// Example: `draft result w 0 0 1 0.2 100`
    pub fn draft(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("draft", args)
    }
    /// DRAW command `drawcont`: 6 uses in 6 converted scripts.
    ///
    /// Example: `drawcont`
    pub fn drawcont(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("drawcont", args)
    }
    /// DRAW command `dset`: 476 uses in 388 converted scripts.
    ///
    /// Example: `dset r sqrt(2)`
    pub fn dset(&mut self, args: &[&str]) -> DrawResult {
        // `dset name value [name value ...]`: numeric variables (values are expressions).
        if args.is_empty() || args.len() % 2 != 0 {
            return Err(DrawError::Failed("usage: dset name value [name value ...]".into()));
        }
        for pair in args.chunks(2) {
            let value = self.number(pair[1])?;
            self.set(pair[0], DrawObject::Number(value));
        }
        Ok(String::new())
    }
    /// DRAW command `dval`: 81 uses in 34 converted scripts.
    ///
    /// Example: `dval dv1_val`
    pub fn dval(&mut self, args: &[&str]) -> DrawResult {
        // `dval expression`: its value.
        need(args, 1, "dval expression")?;
        // (A name no variable has is 0, as in DRAW.)
        let name = args[0].trim();
        if !name.is_empty() && name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') && !name.starts_with(|c: char| c.is_ascii_digit()) && !self.vars.contains_key(name) && name != "pi" {
            return Ok(format_real(0.0));
        }
        Ok(format_real(self.number(args[0])?))
    }
    /// DRAW command `edge`: 875 uses in 170 converted scripts.
    ///
    /// Example: `edge e1 va vb`
    pub fn edge(&mut self, args: &[&str]) -> DrawResult {
        // `edge name v1 v2`: the straight edge between two vertices.
        need(args, 3, "edge name v1 v2")?;
        let (a, b) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let edge = make::edge(&mut self.model, a, b)?;
        self.set_shape_var(args[0], edge)
    }
    /// DRAW command `ellipse`: 96 uses in 92 converted scripts.
    ///
    /// Example: `ellipse w1 0 0 0 15 10`
    pub fn ellipse(&mut self, args: &[&str]) -> DrawResult {
        // `ellipse name x y [dx dy] major minor` (2D) or `ellipse name x y z [nx ny nz [ux uy
        // uz]] major minor`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(2));
        let object = match frame_args.len() {
            2 | 4 => DrawObject::Curve2(Curve2::Ellipse(Ellipse2::new(Draw::curve_frame2(frame_args)?, r[0], r[1]))),
            _ => DrawObject::Curve(Curve::Ellipse(Ellipse::new(Draw::curve_frame(frame_args)?, r[0], r[1]))),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `edgestofaces`.
    ///
    /// Example: `edgestofaces result c`
    pub fn edgestofaces(&mut self, args: &[&str]) -> DrawResult {
        // `edgestofaces result edges [-a angular] [-s shared]`: the planar faces the edges
        // bound (`-s 1`: they already share their vertices; the angular tolerance is not
        // needed).
        need(args, 2, "edgestofaces result edges [-a angle] [-s 0/1]")?;
        let edges = self.shape_var(args[1])?;
        let mut shared = false;
        let mut k = 2;
        while k < args.len() {
            match args[k] {
                "-s" if k + 1 < args.len() => {
                    shared = args[k + 1] != "0";
                    k += 2;
                }
                "-a" if k + 1 < args.len() => k += 2,
                other => return Err(DrawError::Failed(format!("edgestofaces: unknown option {other}"))),
            }
        }
        let result = kentron::make::faces_from_edges(&mut self.model, edges, shared)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `igesbrep`.
    ///
    /// Example: `igesbrep file.igs a *`
    pub fn igesbrep(&mut self, args: &[&str]) -> DrawResult {
        // `igesbrep file name [entities]`: the file's independent entities, a compound named
        // `name` (and kept for `tpcompound`).
        need(args, 2, "igesbrep file name [*]")?;
        let shape = self.read_iges(args[0])?;
        self.last_transfer = Some(shape);
        self.set_shape_var(args[1], shape)
    }
    /// DRAW command `igesread`.
    ///
    /// Example: `igesread file.igs result *`
    pub fn igesread(&mut self, args: &[&str]) -> DrawResult {
        self.igesbrep(args)
    }
    /// DRAW command `testreadiges`.
    ///
    /// Example: `testreadiges file.igs b1`
    pub fn testreadiges(&mut self, args: &[&str]) -> DrawResult {
        need(args, 2, "testreadiges file name")?;
        let shape = self.read_iges(args[0])?;
        self.set_shape_var(args[1], shape)
    }
    /// DRAW command `tpcompound`.
    ///
    /// Example: `tpcompound result`
    pub fn tpcompound(&mut self, args: &[&str]) -> DrawResult {
        // `tpcompound name`: what the last reading transferred, as a compound.
        need(args, 1, "tpcompound name")?;
        let shape = self.last_transfer.ok_or_else(|| DrawError::Failed("tpcompound: nothing transferred".into()))?;
        self.set_shape_var(args[0], shape)
    }
    /// An IGES file (found as `locate_data_file` finds it) read into the model.
    fn read_iges(&mut self, file: &str) -> std::result::Result<kentron::topology::Shape, DrawError> {
        let path = super::find_data_file(file).unwrap_or_else(|| file.to_string());
        let mut f = std::fs::File::open(&path).map_err(|e| DrawError::Failed(format!("{path}: {e}")))?;
        Ok(kentron::io::iges::read(&mut self.model, &mut f)?)
    }
    /// DRAW command `emptycopy`: 5 uses in 5 converted scripts.
    ///
    /// Example: `emptycopy b1 b`
    pub fn emptycopy(&mut self, args: &[&str]) -> DrawResult {
        // `emptycopy result shape`: a shape of the same kind and geometry, without sub-shapes.
        need(args, 2, "emptycopy result shape")?;
        let shape = self.shape_var(args[1])?;
        let copy = self.model.empty_copy(shape);
        self.set_shape_var(args[0], copy)
    }
    /// DRAW command `encoderegularity`: 11 uses in 11 converted scripts.
    ///
    /// Example: `encoderegularity a 1`
    pub fn encoderegularity(&mut self, args: &[&str]) -> DrawResult {
        // `encoderegularity shape [angle]`: records which edges join their faces smoothly. kentron
        // finds that from the geometry where it needs it (as fillets do), so there is nothing to
        // record.
        need(args, 1, "encoderegularity shape [angle]")?;
        self.shape_var(args[0])?;
        Ok(String::new())
    }
    /// DRAW command `evolved`: 8 uses in 8 converted scripts.
    ///
    /// Example: `evolved res -s a -p b -a`
    pub fn evolved(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("evolved", args)
    }
    /// DRAW command `exchuv`: 1 uses in 1 converted scripts.
    ///
    /// Example: `exchuv bs`
    pub fn exchuv(&mut self, args: &[&str]) -> DrawResult {
        // `exchuv surface ...`: the B-spline or Bézier surfaces with u and v exchanged.
        need(args, 1, "exchuv surface ...")?;
        for name in args {
            let s = match self.surface_var(name)? {
                Surface::BSpline(mut b) => { b.exchange_uv(); Surface::BSpline(b) }
                Surface::Bezier(mut b) => { b.exchange_uv(); Surface::Bezier(b) }
                _ => return Err(DrawError::Failed(format!("{name} is not a B-spline or Bézier surface"))),
            };
            self.set(name, DrawObject::Surface(s));
        }
        Ok(String::new())
    }
    /// DRAW command `explode`: 1629 uses in 1315 converted scripts.
    ///
    /// Example: `explode s e`
    pub fn explode(&mut self, args: &[&str]) -> DrawResult {
        // `explode name [C | So | Sh | F | W | E | V]`: the sub-shapes of a type (each once), or
        // the direct children, as variables name_1, name_2, ...
        need(args, 1, "explode name [type]")?;
        let shape = self.shape_var(args[0])?;
        let parts: Vec<Shape> = match args.get(1).map(|t| t.to_ascii_uppercase()) {
            None => self.model.children(shape).collect(),
            Some(t) => {
                let kind = match t.as_str() {
                    "C" | "COMPOUND" => NodeKind::Compound,
                    "CS" | "COMPSOLID" => NodeKind::CompSolid,
                    "SO" | "SOLID" => NodeKind::Solid,
                    "SH" | "SHELL" => NodeKind::Shell,
                    "F" | "FACE" => NodeKind::Face,
                    "W" | "WIRE" => NodeKind::Wire,
                    "E" | "EDGE" => NodeKind::Edge,
                    "V" | "VERTEX" => NodeKind::Vertex,
                    _ => return Err(DrawError::Failed(format!("explode: unknown type {t}"))),
                };
                let mut seen = Vec::new();
                for s in self.model.explore(shape, kind) {
                    if !seen.iter().any(|q: &Shape| q.is_same(&s)) {
                        seen.push(s);
                    }
                }
                seen
            }
        };
        let names: Vec<String> = (1..=parts.len()).map(|i| format!("{}_{i}", args[0])).collect();
        for (name, part) in names.iter().zip(parts) {
            self.set_shape_var(name, part)?;
        }
        Ok(names.join(" "))
    }
    /// DRAW command `extrema`: 9 uses in 8 converted scripts.
    ///
    /// Example: `extrema result l`
    pub fn extrema(&mut self, args: &[&str]) -> DrawResult {
        // `extrema curve|surface curve|surface`: the extremal distances between a curve and a
        // curve or surface, each as `ext_i`, the segment between the extremal points (none where
        // they coincide); the distance alone where the two run parallel.
        need(args, 2, "extrema curve|surface curve|surface")?;
        use kentron::geometry::extrema::{CurveCurveExtrema, CurveSurfaceExtrema};
        let (a, b) = (self.get(args[0]).cloned(), self.get(args[1]).cloned());
        let pairs: Vec<(Point, Point)> = match (a, b) {
            (Some(DrawObject::Curve(c1)), Some(DrawObject::Curve(c2))) => match kentron::geometry::extrema::curve_curve(&c1, &c2) {
                CurveCurveExtrema::Isolated(v) => v.into_iter().map(|e| (e.p1, e.p2)).collect(),
                CurveCurveExtrema::Parallel { distance } => return Ok(format!("Infinite number of extremas, distance = {}", format_real(distance))),
            },
            (Some(DrawObject::Curve(c)), Some(DrawObject::Surface(s))) => match kentron::geometry::extrema::curve_surface(&c, &s) {
                CurveSurfaceExtrema::Isolated(v) => v.into_iter().map(|e| (e.curve_point, e.surface_point)).collect(),
                CurveSurfaceExtrema::Parallel { distance } => return Ok(format!("Infinite number of extremas, distance = {}", format_real(distance))),
            },
            (Some(DrawObject::Surface(s)), Some(DrawObject::Curve(c))) => match kentron::geometry::extrema::curve_surface(&c, &s) {
                CurveSurfaceExtrema::Isolated(v) => v.into_iter().map(|e| (e.surface_point, e.curve_point)).collect(),
                CurveSurfaceExtrema::Parallel { distance } => return Ok(format!("Infinite number of extremas, distance = {}", format_real(distance))),
            },
            _ => return self.not_implemented("extrema", args),
        };
        if pairs.is_empty() {
            return Ok("No solutions!".into());
        }
        let mut out = String::new();
        for (i, (p, q)) in pairs.into_iter().enumerate() {
            let d = p.distance(q);
            if d < 1e-16 {
                out += &format!("Extrema {} is point : {} {} {}\n", i + 1, p.x, p.y, p.z);
                continue;
            }
            let line = Line::new(p, Direction::from_vector(q - p));
            let name = format!("ext_{}", i + 1);
            self.set(&name, DrawObject::Curve(Curve::Trimmed(TrimmedCurve::new(Curve::Line(line), 0.0, d))));
            out += &format!("{name} ");
        }
        Ok(out)
    }
    /// DRAW command `extsurf`: 52 uses in 52 converted scripts.
    ///
    /// Example: `extsurf ss cc 1 0 0`
    pub fn extsurf(&mut self, args: &[&str]) -> DrawResult {
        // `extsurf name curve dx dy dz`: the surface swept by the curve (or an edge's) along
        // a direction.
        need(args, 5, "extsurf name curve dx dy dz")?;
        let curve = self.curve_or_edge(args[1])?;
        let v = self.numbers(&args[2..5])?;
        let d = Direction::try_from_vector(Vector::new(v[0], v[1], v[2])).ok_or_else(|| DrawError::Failed("extsurf: null direction".into()))?;
        self.set(args[0], DrawObject::Surface(Surface::Extrusion(ExtrusionSurface::new(curve, d))));
        Ok(String::new())
    }
    /// DRAW command `featdprism`: 36 uses in 36 converted scripts.
    ///
    /// Example: `featdprism s a a 25 0 1`
    pub fn featdprism(&mut self, args: &[&str]) -> DrawResult {
        // `featdprism base profile sketch-face angle(degrees) fuse modify`.
        let [base, profile, _sketch, angle, fuse, _modify] = args else { return self.not_implemented("featdprism", args) };
        // (DRAW converts the angle with π/360, halving it: kept, to match its results.)
        let angle = self.number(angle)? * std::f64::consts::PI / 360.0;
        self.define_feature("dprism", base, profile, kentron::feature::FeatureSweep::DraftedPrism(angle), fuse)
    }
    /// DRAW command `featlf`: 41 uses in 40 converted scripts.
    ///
    /// Example: `featlf B w p -5 0 0 5 0 0 1 1`
    pub fn featlf(&mut self, args: &[&str]) -> DrawResult {
        // `featlf base wire plane dx dy dz dx1 dy1 dz1 fuse modify`: the wire thickened to the
        // first vector on one side and the second on the other.
        if args.len() != 11 {
            return self.not_implemented("featlf", args);
        }
        let v = self.numbers(&args[3..9])?;
        let sweep = kentron::feature::RibSweep::Linear { side: Vector::new(v[0], v[1], v[2]), other_side: Vector::new(v[3], v[4], v[5]) };
        self.define_rib("lf", args[0], args[1], args[2], sweep, args[9])
    }
    /// DRAW command `featperform`: 220 uses in 218 converted scripts.
    ///
    /// Example: `featperform revol result boite_1`
    pub fn featperform(&mut self, args: &[&str]) -> DrawResult {
        // `featperform kind result [until | from until]`: through all, up to a face (to the
        // end of the shape if it is not a shape), or between two faces.
        use kentron::feature::FeatureExtent;
        let (kind, result) = match args {
            [k, r, ..] => (*k, *r),
            _ => return self.not_implemented("featperform", args),
        };
        let face = |d: &Self, name: &str| match d.get(name) { Some(DrawObject::Shape(s)) => Some(*s), _ => None };
        let extent = match &args[2..] {
            [] => FeatureExtent::ThruAll,
            [until] => match face(self, until) { Some(u) => FeatureExtent::Until(u), None => FeatureExtent::UntilEnd },
            [from, until] => match (face(self, from), face(self, until)) {
                (Some(f), Some(u)) => FeatureExtent::FromUntil(f, u),
                _ => return self.not_implemented("featperform", args),
            },
            _ => return self.not_implemented("featperform", args),
        };
        self.perform_feature(kind, result, extent)
    }
    /// DRAW command `featperformval`: 93 uses in 92 converted scripts.
    ///
    /// Example: `featperformval dprism result 200`
    pub fn featperformval(&mut self, args: &[&str]) -> DrawResult {
        // `featperformval kind result value [until]`: a length (an angle in degrees for a
        // revolution), towards a face and stopping there if given.
        let (kind, result, value, until) = match args {
            [k, r, v] => (*k, *r, *v, None),
            [k, r, v, u] => (*k, *r, *v, Some(self.shape_var(u)?)),
            _ => return self.not_implemented("featperformval", args),
        };
        let mut value = self.number(value)?;
        if kind.eq_ignore_ascii_case("revol") {
            value = value.to_radians();
        }
        let extent = match until {
            Some(u) => kentron::feature::FeatureExtent::UntilLength(u, value),
            None => kentron::feature::FeatureExtent::Length(value),
        };
        self.perform_feature(kind, result, extent)
    }
    /// Keeps a feature's definition for `featperform`.
    fn define_feature(&mut self, kind: &str, base: &str, profile: &str, sweep: kentron::feature::FeatureSweep, fuse: &str) -> DrawResult {
        use kentron::feature::FeatureMode;
        let base = self.shape_var(base)?;
        let profile = self.shape_var(profile)?;
        let mode = match self.number(fuse)? as i64 {
            0 => FeatureMode::Remove,
            1 => FeatureMode::Add,
            _ => FeatureMode::ToolOnly,
        };
        self.features.insert(kind.to_ascii_lowercase(), (base, profile, sweep, mode));
        Ok(String::new())
    }

    /// Keeps a rib's or slot's definition for `featperform`.
    fn define_rib(&mut self, kind: &str, base: &str, wire: &str, plane: &str, sweep: kentron::feature::RibSweep, fuse: &str) -> DrawResult {
        use kentron::feature::FeatureMode;
        let base = self.shape_var(base)?;
        let wire = self.shape_var(wire)?;
        let plane = match self.get(plane) {
            Some(DrawObject::Surface(Surface::Plane(p))) => *p,
            _ => return Err(DrawError::Failed(format!("{plane} is not a plane"))),
        };
        let mode = match self.number(fuse)? as i64 {
            0 => FeatureMode::Remove,
            1 => FeatureMode::Add,
            _ => FeatureMode::ToolOnly,
        };
        self.ribs.insert(kind.to_string(), (base, wire, plane, sweep, mode));
        Ok(String::new())
    }

    /// Makes the feature of kind `kind` defined last, to `extent`, as `result`.
    /// (Run on a copy of the model in a thread given 60 s, like the YAML cases: some booleans
    /// on filleted bases do not finish.)
    fn perform_feature(&mut self, kind: &str, result: &str, extent: kentron::feature::FeatureExtent) -> DrawResult {
        let rib = self.ribs.get(&kind.to_ascii_lowercase()).copied();
        let feature = self.features.get(&kind.to_ascii_lowercase()).copied();
        if rib.is_none() && feature.is_none() {
            return Err(DrawError::Failed(format!("featperform: no {kind} defined")));
        }
        let mut model = self.model.clone();
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let shape = match (rib, feature) {
                (Some((base, wire, plane, sweep, mode)), _) => kentron::feature::rib(&mut model, base, wire, &plane, sweep, mode),
                (_, Some((base, profile, sweep, mode))) => kentron::feature::feature(&mut model, base, profile, sweep, mode, extent),
                _ => unreachable!(),
            };
            let _ = tx.send((shape, model));
        });
        match rx.recv_timeout(std::time::Duration::from_secs(60)) {
            Ok((shape, model)) => {
                self.model = model;
                self.set_shape_var(result, shape?)
            }
            Err(std::sync::mpsc::RecvTimeoutError::Timeout) => Err(DrawError::Failed("featperform: timed out after 60 s".into())),
            Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => panic!("featperform panicked"),
        }
    }
    /// DRAW command `featprism`: 143 uses in 141 converted scripts.
    ///
    /// Example: `featprism aB aB aB 0 -1 0 2 1`
    pub fn featprism(&mut self, args: &[&str]) -> DrawResult {
        // `featprism base profile sketch-face dx dy dz fuse(0 remove/1 add/2 tool only) modify`.
        let [base, profile, _sketch, dx, dy, dz, fuse, _modify] = args else { return self.not_implemented("featprism", args) };
        let v = self.numbers(&[dx, dy, dz])?;
        let direction = Direction::try_from_vector(Vector::new(v[0], v[1], v[2])).ok_or_else(|| DrawError::Failed("featprism: null direction".into()))?;
        self.define_feature("prism", base, profile, kentron::feature::FeatureSweep::Prism(direction), fuse)
    }
    /// DRAW command `featrevol`: 80 uses in 80 converted scripts.
    ///
    /// Example: `featrevol boite cont boite_4 0 0 200 0 -1 0 1 1`
    pub fn featrevol(&mut self, args: &[&str]) -> DrawResult {
        // `featrevol base profile sketch-face ox oy oz dx dy dz fuse modify`.
        if args.len() != 11 {
            return self.not_implemented("featrevol", args);
        }
        let v = self.numbers(&args[3..9])?;
        let direction = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("featrevol: null axis".into()))?;
        let axis = Axis::new(Point::new(v[0], v[1], v[2]), direction);
        self.define_feature("revol", args[0], args[1], kentron::feature::FeatureSweep::Revolution(axis), args[9])
    }
    /// DRAW command `featrf`: 13 uses in 13 converted scripts.
    ///
    /// Example: `featrf s w pl 0 0 0 0 0 1 0.2 0.2 1 1`
    pub fn featrf(&mut self, args: &[&str]) -> DrawResult {
        // `featrf base wire plane ox oy oz dx dy dz size1 size2 fuse modify`: the wire turned
        // about the axis, as thick as the sizes on each side of its plane.
        if args.len() != 13 {
            return self.not_implemented("featrf", args);
        }
        let v = self.numbers(&args[3..11])?;
        let direction = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("featrf: null axis".into()))?;
        let sweep = kentron::feature::RibSweep::Revolution { axis: Axis::new(Point::new(v[0], v[1], v[2]), direction), side: v[6], other_side: v[7] };
        self.define_rib("rf", args[0], args[1], args[2], sweep, args[11])
    }
    /// DRAW command `fillcurves`: 2 uses in 1 converted scripts.
    ///
    /// Example: `fillcurves s11 c1 c41 cv c31`
    pub fn fillcurves(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("fillcurves", args)
    }
    /// DRAW command `fillet`: 7 uses in 7 converted scripts.
    ///
    /// Example: `fillet result b1 5 b1_15 5 b1_16 5 b1_10`
    pub fn fillet(&mut self, args: &[&str]) -> DrawResult {
        // `fillet result shape radius edge [radius edge ...]`, or one radius for several edges
        // (`fillet result shape radius e1 e2 ...` is OCCT's older form, one radius per edge
        // following it until the next number).
        need(args, 4, "fillet result shape radius edge [radius edge ...]")?;
        let shape = self.shape_var(args[1])?;
        let mut edges = Vec::new();
        let mut radius = None;
        for a in &args[2..] {
            match self.get(a) {
                Some(DrawObject::Shape(e)) => {
                    let r = radius.ok_or_else(|| DrawError::Failed("fillet: an edge before any radius".into()))?;
                    for edge in self.model.unique_subshapes(*e, NodeKind::Edge) {
                        edges.push((edge, r));
                    }
                }
                _ => radius = Some(self.number(a)?),
            }
        }
        let result = blend::fillet_edges(&mut self.model, shape, &edges)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `fillet2d`: 9 uses in 9 converted scripts.
    ///
    /// Example: `fillet2d result e1 e2 10`
    pub fn fillet2d(&mut self, args: &[&str]) -> DrawResult {
        // `fillet2d result edge1 edge2 radius`: the fillet arc (the edges cut back to it), or
        // `fillet2d result wire radius`: the wire of two edges with its corner rounded.
        self.corner_2d("fillet2d", args, 1)
    }
    /// DRAW command `firsthole`: 2 uses in 2 converted scripts.
    ///
    /// Example: `firsthole r a 7.31762745781211 0 22.6223587092621 -0.309016994374947 -0 -0.951056516295154 5`
    pub fn firsthole(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("firsthole", args)
    }
    /// DRAW command `fixshape`: 87 uses in 74 converted scripts.
    ///
    /// Example: `fixshape b b`
    pub fn fixshape(&mut self, args: &[&str]) -> DrawResult {
        // `fixshape result shape [precision [max-tolerance]] [-maxtaila angle -maxtailw
        // width]`: with the tail options, the wires' tails (narrower than the width, at an
        // angle below the given degrees) are cut off first.
        let mut args: Vec<&str> = args.to_vec();
        let mut tails = (None, None);
        while let Some(k) = args.iter().position(|a| *a == "-maxtaila" || *a == "-maxtailw") {
            let value = self.number(args.get(k + 1).copied().unwrap_or(""))?;
            if args[k] == "-maxtaila" { tails.0 = Some(value) } else { tails.1 = Some(value) }
            args.drain(k..k + 2);
        }
        // (`+o`, `-o` and the like switch fixing modes on and off: kentron fixes all it can.)
        let flag = |a: &str| (a.starts_with('+') || a.starts_with('-')) && a.parse::<f64>().is_err() && a.len() > 1;
        let args: Vec<&str> = args.iter().copied().filter(|a| !flag(a)).collect();
        let args = &args[..];
        if args.len() < 2 || args.len() > 4 {
            return self.not_implemented("fixshape", args);
        }
        let shape = self.shape_var(args[1])?;
        let shape = match tails {
            (None, None) => shape,
            (angle, Some(width)) => {
                let copy = modify::copy(&mut self.model, shape, false)?;
                kentron::heal::fix_wire_tails(&mut self.model, copy, angle.unwrap_or(0.0).to_radians(), width)?
            }
            (Some(_), None) => return self.not_implemented("fixshape", args),
        };
        let mut options = kentron::heal::FixOptions::default();
        if let Some(p) = args.get(2) {
            options.precision = self.number(p)?;
        }
        if let Some(m) = args.get(3) {
            options.max_tolerance = self.number(m)?;
        }
        let result = kentron::heal::fix_shape(&mut self.model, shape, &options)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `fixsmall`: 1 uses in 1 converted scripts.
    ///
    /// Example: `fixsmall result a 0.05`
    pub fn fixsmall(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("fixsmall", args)
    }
    /// DRAW command `fixsmalledges`: 1 uses in 1 converted scripts.
    ///
    /// Example: `fixsmalledges a a 0.002`
    pub fn fixsmalledges(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("fixsmalledges", args)
    }
    /// DRAW command `fixsmallfaces`: 3 uses in 2 converted scripts.
    ///
    /// Example: `fixsmallfaces result shape 1.e-7`
    pub fn fixsmallfaces(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("fixsmallfaces", args)
    }
    /// DRAW command `fsameparameter`: 55 uses in 49 converted scripts.
    ///
    /// Example: `fsameparameter s SCALE1*1.e-7`
    pub fn fsameparameter(&mut self, args: &[&str]) -> DrawResult {
        // `fsameparameter shape [tolerance]`: makes pcurves agree with 3D curves. Kentron's
        // pcurves are built to agree, so there is nothing to do.
        need(args, 1, "fsameparameter shape [tolerance]")?;
        self.shape_var(args[0])?;
        Ok(String::new())
    }
    /// DRAW command `fuseedge`: 1 uses in 1 converted scripts.
    ///
    /// Example: `fuseedge result`
    pub fn fuseedge(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("fuseedge", args)
    }
    /// DRAW command `generated`: 445 uses in 42 converted scripts.
    ///
    /// Example: `generated r1 hist w1_1`
    pub fn generated(&mut self, args: &[&str]) -> DrawResult {
        // `generated result history shape`: a compound of the shapes generated from `shape`.
        need(args, 3, "generated result history shape")?;
        let generated = match self.get(args[1]) {
            Some(DrawObject::History(h)) => h.generated.clone(),
            _ => return Err(DrawError::Failed(format!("{} is not a history", args[1]))),
        };
        let mut input = self.shape_var(args[2])?;
        // (A wire's child, from `explode`, is a coedge: its edge.)
        if input.kind() == NodeKind::CoEdge {
            input = self.model.explore(input, NodeKind::Edge).next().unwrap_or(input);
        }
        let shapes: Vec<Shape> = generated.iter().filter(|(node, _)| *node == input.node).flat_map(|(_, s)| s.clone()).collect();
        if shapes.is_empty() {
            return Ok("No shapes were generated from the shape\n".into());
        }
        let compound = self.compound_of(&shapes)?;
        self.set_shape_var(args[0], compound)
    }
    /// DRAW command `gplate`: 1 uses in 1 converted scripts.
    ///
    /// Example: `gplate result 5 0 c_1 0 c_2 0 c_3 0 c_4 0 c_5 0`
    pub fn gplate(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("gplate", args)
    }
    /// DRAW command `halfspace`: 26 uses in 23 converted scripts.
    ///
    /// Example: `halfspace space2 bsurf2 0 -0.3 0`
    pub fn halfspace(&mut self, args: &[&str]) -> DrawResult {
        // `halfspace result face x y z`: the half-space bounded by the face (or shell) on the
        // side of the point.
        need(args, 5, "halfspace result face x y z")?;
        let face = self.shape_var(args[1])?;
        let v = self.numbers(&args[2..5])?;
        let result = kentron::primitives::half_space(&mut self.model, face, Point::new(v[0], v[1], v[2]))?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `hfill`: 5 uses in 4 converted scripts.
    ///
    /// Example: `hfill Shape2DView_outl Shape2DView_proj 0`
    pub fn hfill(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hfill", args)
    }
    /// DRAW command `hhide`: 4 uses in 4 converted scripts.
    ///
    /// Example: `hhide`
    pub fn hhide(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hhide", args)
    }
    /// DRAW command `hload`: 5 uses in 4 converted scripts.
    ///
    /// Example: `hload Shape2DView_outl`
    pub fn hload(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hload", args)
    }
    /// DRAW command `hlr`: 2 uses in 1 converted scripts.
    ///
    /// Example: `hlr hlr a`
    pub fn hlr(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hlr", args)
    }
    /// DRAW command `hlrin3d`: 1 uses in 1 converted scripts.
    ///
    /// Example: `hlrin3d result a 0.57735026918962573 -0.57735026918962573 0.57735026918962584`
    pub fn hlrin3d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hlrin3d", args)
    }
    /// DRAW command `houtl`: 5 uses in 4 converted scripts.
    ///
    /// Example: `houtl Shape2DView_outl Line`
    pub fn houtl(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("houtl", args)
    }
    /// DRAW command `hprj`: 4 uses in 4 converted scripts.
    ///
    /// Example: `hprj Shape2DView_proj 0 0 0 0 0 1 1 0 0`
    pub fn hprj(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hprj", args)
    }
    /// DRAW command `hremove`: 1 uses in 1 converted scripts.
    ///
    /// Example: `hremove`
    pub fn hremove(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hremove", args)
    }
    /// DRAW command `hres2d`: 4 uses in 4 converted scripts.
    ///
    /// Example: `hres2d`
    pub fn hres2d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hres2d", args)
    }
    /// DRAW command `hsetprj`: 4 uses in 4 converted scripts.
    ///
    /// Example: `hsetprj Shape2DView_proj`
    pub fn hsetprj(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hsetprj", args)
    }
    /// DRAW command `hupdate`: 4 uses in 4 converted scripts.
    ///
    /// Example: `hupdate`
    pub fn hupdate(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("hupdate", args)
    }
    /// DRAW command `hyperbola`: 41 uses in 40 converted scripts.
    ///
    /// Example: `hyperbola result 0 0 0 0 0 1 1 0 0 20 10`
    pub fn hyperbola(&mut self, args: &[&str]) -> DrawResult {
        // `hyperbola name x y [dx dy] major minor` (2D) or `hyperbola name x y z [nx ny nz [ux
        // uy uz]] major minor`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(2));
        let object = match frame_args.len() {
            2 | 4 => DrawObject::Curve2(Curve2::Hyperbola(Hyperbola2::new(Draw::curve_frame2(frame_args)?, r[0], r[1]))),
            _ => DrawObject::Curve(Curve::Hyperbola(Hyperbola::new(Draw::curve_frame(frame_args)?, r[0], r[1]))),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `incdeg`: 5 uses in 5 converted scripts.
    ///
    /// Example: `incdeg result 8`
    pub fn incdeg(&mut self, args: &[&str]) -> DrawResult {
        // `incdeg curve degree`: the B-spline or Bézier curve (3D or 2D) raised to the degree.
        need(args, 2, "incdeg curve degree")?;
        let degree = self.number(args[1])? as usize;
        let raised = match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(Curve::BSpline(mut b))) => DrawObject::Curve(Curve::BSpline({ b.increase_degree(degree); b })),
            Some(DrawObject::Curve(Curve::Bezier(mut b))) => DrawObject::Curve(Curve::Bezier({ b.increase_degree(degree); b })),
            Some(DrawObject::Curve2(Curve2::BSpline(mut b))) => DrawObject::Curve2(Curve2::BSpline({ b.increase_degree(degree); b })),
            Some(DrawObject::Curve2(Curve2::Bezier(mut b))) => DrawObject::Curve2(Curve2::Bezier({ b.increase_degree(degree); b })),
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline or Bézier curve", args[0]))),
        };
        self.set(args[0], raised);
        Ok(String::new())
    }
    /// DRAW command `incmesh`: 169 uses in 153 converted scripts.
    ///
    /// Example: `incmesh result 0.01`
    pub fn incmesh(&mut self, args: &[&str]) -> DrawResult {
        // `incmesh shape... deflection [-a angle_degrees] [other options]`: meshes the shapes
        // (the first number is the deflection; `-relative` and the other options are not
        // taken into account).
        let mut shapes = Vec::new();
        let mut deflection = None;
        let mut angle = 0.5;
        let mut i = 0;
        while i < args.len() {
            match args[i] {
                "-a" => {
                    angle = self.number(args.get(i + 1).copied().unwrap_or("0.5"))?.to_radians();
                    i += 1;
                }
                // (Options with a value, not taken into account.)
                "-ai" | "-min" | "-algo" | "-di" => i += 1,
                a if a.starts_with('-') => {}
                a if deflection.is_none() && self.number(a).is_ok() && self.shape_var(a).is_err() => deflection = Some(self.number(a)?),
                // (A number after the deflection, as the value OCCT's scripts give `-parallel`:
                // not taken into account.)
                a if deflection.is_some() && self.number(a).is_ok() && self.shape_var(a).is_err() => {}
                a => shapes.push(self.shape_var(a)?),
            }
            i += 1;
        }
        let deflection = deflection.ok_or_else(|| DrawError::Failed("incmesh: no deflection".into()))?;
        for shape in shapes {
            kentron::meshing::mesh_with(&mut self.model, shape, kentron::meshing::MeshOptions { deflection, angle })?;
        }
        Ok(String::new())
    }
    /// DRAW command `incudeg`: 4 uses in 4 converted scripts.
    ///
    /// Example: `incudeg p 3`
    pub fn incudeg(&mut self, args: &[&str]) -> DrawResult {
        // `incudeg surface degree`: the B-spline or Bézier surface raised to the degree in u.
        need(args, 2, "incudeg surface degree")?;
        let degree = self.number(args[1])? as usize;
        let raised = match self.surface_var(args[0])? {
            Surface::BSpline(mut b) => { let v = b.v_degree(); b.increase_degree(degree, v); Surface::BSpline(b) }
            Surface::Bezier(mut b) => { let v = b.v_degree(); b.increase_degree(degree, v); Surface::Bezier(b) }
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline or Bézier surface", args[0]))),
        };
        self.set(args[0], DrawObject::Surface(raised));
        Ok(String::new())
    }
    /// DRAW command `incvdeg`: 4 uses in 4 converted scripts.
    ///
    /// Example: `incvdeg p 3`
    pub fn incvdeg(&mut self, args: &[&str]) -> DrawResult {
        // `incvdeg surface degree`: the B-spline or Bézier surface raised to the degree in v.
        need(args, 2, "incvdeg surface degree")?;
        let degree = self.number(args[1])? as usize;
        let raised = match self.surface_var(args[0])? {
            Surface::BSpline(mut b) => { let u = b.u_degree(); b.increase_degree(u, degree); Surface::BSpline(b) }
            Surface::Bezier(mut b) => { let u = b.u_degree(); b.increase_degree(u, degree); Surface::Bezier(b) }
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline or Bézier surface", args[0]))),
        };
        self.set(args[0], DrawObject::Surface(raised));
        Ok(String::new())
    }
    /// DRAW command `insertknot`: 2 uses in 2 converted scripts.
    ///
    /// Example: `insertknot result 0.5 1`
    pub fn insertknot(&mut self, args: &[&str]) -> DrawResult {
        // `insertknot curve knot [mult] [knot mult ...]`: knots inserted into the B-spline curve
        // (3D or 2D), each to the multiplicity given (1 by default).
        need(args, 2, "insertknot curve knot [mult] ...")?;
        let v = self.numbers(&args[1..])?;
        let pairs: Vec<(f64, usize)> = if v.len() == 1 { vec![(v[0], 1)] } else { v.chunks(2).map(|c| (c[0], c.get(1).copied().unwrap_or(1.0) as usize)).collect() };
        match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(Curve::BSpline(mut b))) => {
                for &(u, m) in &pairs { b.insert_knot(u, m); }
                self.set(args[0], DrawObject::Curve(Curve::BSpline(b)));
            }
            Some(DrawObject::Curve2(Curve2::BSpline(mut b))) => {
                for &(u, m) in &pairs { b.insert_knot(u, m); }
                self.set(args[0], DrawObject::Curve2(Curve2::BSpline(b)));
            }
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline curve", args[0]))),
        }
        Ok(String::new())
    }
    /// DRAW command `insertpole`: 2 uses in 2 converted scripts.
    ///
    /// Example: `insertpole result 2 1 6 1`
    pub fn insertpole(&mut self, args: &[&str]) -> DrawResult {
        // `insertpole curve index x y [z] [weight]`: a pole inserted after the index (from 1)
        // into the Bézier curve (2D: x y [weight]).
        need(args, 4, "insertpole curve index x y [z] [weight]")?;
        // (DRAW counts poles from 1: after pole `index`.)
        let index = (self.number(args[1])? as usize).saturating_sub(1);
        let v = self.numbers(&args[2..])?;
        match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(Curve::Bezier(mut b))) if v.len() == 3 || v.len() == 4 => {
                let p = Point::new(v[0], v[1], v[2]);
                match v.get(3) {
                    Some(&w) => b.insert_pole_after_with_weight(index, p, w),
                    None => b.insert_pole_after(index, p),
                }
                self.set(args[0], DrawObject::Curve(Curve::Bezier(b)));
            }
            Some(DrawObject::Curve2(Curve2::Bezier(mut b))) if v.len() == 2 || v.len() == 3 => {
                let p = Point2::new(v[0], v[1]);
                match v.get(2) {
                    Some(&w) => b.insert_pole_after_with_weight(index, p, w),
                    None => b.insert_pole_after(index, p),
                }
                self.set(args[0], DrawObject::Curve2(Curve2::Bezier(b)));
            }
            _ => return self.not_implemented("insertpole", args),
        }
        Ok(String::new())
    }
    /// DRAW command `intersect`: 2 uses in 2 converted scripts.
    ///
    /// Example: `intersect result s1 s2`
    pub fn intersect(&mut self, args: &[&str]) -> DrawResult {
        // `intersect name curve|surface surface [u1 v1 u2 v2] [bounds of both] [tolerance]`:
        // the curves and points where they meet. Two surfaces' curves are named `name` if
        // one, else `name_1`, `name_2`, ..., their points `name_p_1`, ...; a curve's segments
        // and points on a surface `name` if one, else `name_1`, ... Returns the names.
        use kentron::geometry::intersect as isect;
        need(args, 3, "intersect name curve|surface surface [u1 v1 u2 v2] [bounds] [tolerance]")?;
        let n = args.len();
        if !matches!(n, 3 | 4 | 7 | 8 | 11 | 12 | 15 | 16) {
            return Err(DrawError::Failed("intersect: incorrect number of arguments".into()));
        }
        let tolerance = if n % 2 == 0 { self.number(args[n - 1])? } else { 1e-7 };
        let name = args[0];
        let numbered = |count: usize, k: usize| if count == 1 { name.to_string() } else { format!("{name}_{}", k + 1) };
        if let Some(DrawObject::Curve(curve)) = self.get(args[1]).cloned() {
            let surface = self.surface_var(args[2])?;
            let found = isect::curve_surface(&curve, &surface, tolerance);
            let mut objects: Vec<DrawObject> = found.segments.iter().map(|&(a, b)| DrawObject::Curve(Curve::Trimmed(TrimmedCurve::new(curve.clone(), a, b)))).collect();
            objects.extend(found.points.iter().map(|p| DrawObject::Point(p.point)));
            let mut names = Vec::new();
            for (k, object) in objects.iter().enumerate() {
                names.push(numbered(objects.len(), k));
                self.set(names.last().unwrap(), object.clone());
            }
            return Ok(names.join(" "));
        }
        let (s1, s2) = (self.surface_var(args[1])?, self.surface_var(args[2])?);
        let bounds = match n {
            11 | 12 => Some(self.numbers(&args[3..11])?),
            15 | 16 => Some(self.numbers(&args[7..15])?),
            _ => None,
        };
        let (curves, points) = match bounds {
            None => match elementary_intersection(&s1, &s2, tolerance) {
                Some(found) => found,
                None => (self.intersect_by_section(&s1, &s2, None)?, Vec::new()),
            },
            Some(b) => (self.intersect_by_section(&s1, &s2, Some(&b))?, Vec::new()),
        };
        let mut names = Vec::new();
        for (k, curve) in curves.iter().enumerate() {
            names.push(numbered(curves.len(), k));
            self.set(names.last().unwrap(), DrawObject::Curve(curve.clone()));
        }
        for (k, p) in points.iter().enumerate() {
            names.push(format!("{name}_p_{}", k + 1));
            self.set(names.last().unwrap(), DrawObject::Point(*p));
        }
        Ok(names.join(" "))
    }

    /// The curves where two surfaces meet, as a section of faces on their domains (or on
    /// `bounds`, both surfaces' `u1 u2 v1 v2`), unbounded ones limited to a large box;
    /// fails if a curve runs out to it (OCCT finds no intersection then).
    fn intersect_by_section(&mut self, s1: &Surface, s2: &Surface, bounds: Option<&[f64]>) -> DrawResult<Vec<Curve>> {
        const FAR: f64 = 1e4;
        let mut faces = Vec::new();
        for (k, s) in [s1, s2].into_iter().enumerate() {
            let (mut u0, mut u1, mut v0, mut v1) = kentron::geometry::ParametricSurface::bounds(s);
            if let Some(b) = bounds {
                (u0, u1, v0, v1) = (b[4 * k], b[4 * k + 1], b[4 * k + 2], b[4 * k + 3]);
            }
            let clamp = |a: f64, b: f64| if a.abs() > 1e100 || b.abs() > 1e100 { (a.max(-FAR), b.min(FAR)) } else { (a, b) };
            ((u0, u1), (v0, v1)) = (clamp(u0, u1), clamp(v0, v1));
            faces.push(make::face_on_domain(&mut self.model, s.clone(), u0, u1, v0, v1)?);
        }
        let section = boolean::section(&mut self.model, faces[0], faces[1])?;
        let edges = self.model.unique_subshapes(section, NodeKind::Edge);
        for &e in &edges {
            if let Some((c, first, last)) = self.model.edge_curve(e) {
                if [first, last].iter().any(|&t| {
                    let p = c.point(t);
                    [p.x, p.y, p.z].iter().any(|x: &f64| x.abs() > FAR / 2.0)
                }) {
                    return Err(DrawError::Failed("intersect: an intersection curve is unbounded".into()));
                }
            }
        }
        chain_edges(&self.model, &edges).iter().map(|chain| self.fit_on_surfaces(chain, s1, s2)).collect()
    }

    /// A curve through a chain of section edges, its points moved exactly onto both
    /// surfaces, and more of them where it bends.
    fn fit_on_surfaces(&self, chain: &[Shape], s1: &Surface, s2: &Surface) -> DrawResult<Curve> {
        let mut points: Vec<Point> = Vec::new();
        for &e in chain {
            let Some((c, first, last)) = self.model.edge_curve(e) else { continue };
            let (from, to) = if e.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
            let n = 32;
            for i in 0..=n {
                points.push(c.point(from + (to - from) * i as f64 / n as f64));
            }
        }
        match fit_through_on_surfaces(&points, s1, s2) {
            Some(curve) => Ok(curve),
            None => chain_curve(&self.model, chain),
        }
    }
    /// DRAW command `invert`: 8 uses in 8 converted scripts.
    ///
    /// Example: `invert b`
    pub fn invert(&mut self, args: &[&str]) -> DrawResult {
        // `invert shape`: the shape with its orientation reversed.
        need(args, 1, "invert shape")?;
        let shape = self.shape_var(args[0])?;
        self.set_shape_var(args[0], shape.reversed())
    }
    /// DRAW command `l`: 7 uses in 7 converted scripts.
    ///
    /// Example: `l`
    pub fn l(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("l", args)
    }
    /// DRAW command `law`: 4 uses in 4 converted scripts.
    ///
    /// Example: `law law1 5 4 0.0e0 6 1.0e0 1 2.0e0 1 3.0e0 6 0 0 0.1 0.2 0.7 0.3 1 1`
    pub fn law(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("law", args)
    }
    /// DRAW command `line`: 85 uses in 43 converted scripts.
    ///
    /// Example: `line l 0 0 0 0 1 0`
    pub fn line(&mut self, args: &[&str]) -> DrawResult {
        // `line name x y z dx dy dz` or `line name x y dx dy` (2D).
        need(args, 5, "line name x y [z] dx dy [dz]")?;
        let v = self.numbers(&args[1..])?;
        let object = match v.len() {
            4 => {
                let d = Direction2::try_new(v[2], v[3]).ok_or_else(|| DrawError::Failed("line: null direction".into()))?;
                DrawObject::Curve2(Curve2::Line(Line2::new(Point2::new(v[0], v[1]), d)))
            }
            6 => {
                let d = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("line: null direction".into()))?;
                DrawObject::Curve(Curve::Line(Line::new(Point::new(v[0], v[1], v[2]), d)))
            }
            _ => return Err(DrawError::Failed("line: syntax error".into())),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `lmirror`: 9 uses in 8 converted scripts.
    ///
    /// Example: `lmirror result 0 0 0 1 0 0`
    pub fn lmirror(&mut self, args: &[&str]) -> DrawResult {
        // `lmirror name ... x y z dx dy dz`: the symmetry about an axis.
        need(args, 7, "lmirror name ... x y z dx dy dz")?;
        let n = args.len();
        let v = self.numbers(&args[n - 6..])?;
        let d = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("lmirror: null direction".into()))?;
        self.transform_vars(&args[..n - 6], &Transform::mirror_axis(Axis::new(Point::new(v[0], v[1], v[2]), d)))
    }
    /// DRAW command `sprops`: the area of a shape's faces and their centre of mass.
    ///
    /// Example: `sprops result`
    pub fn sprops(&mut self, args: &[&str]) -> DrawResult {
        // `sprops shape [eps]`: the area of its faces and their centre of mass.
        need(args, 1, "sprops shape [eps]")?;
        let shape = self.shape_var(args[0])?;
        Ok(mass_report(properties::area(&self.model, shape)))
    }
    /// DRAW command `vprops`: the volume of a shape and its centre of mass.
    ///
    /// Example: `vprops s`
    pub fn vprops(&mut self, args: &[&str]) -> DrawResult {
        // `vprops shape [eps]`
        need(args, 1, "vprops shape [eps]")?;
        let shape = self.shape_var(args[0])?;
        Ok(mass_report(properties::volume(&self.model, shape)))
    }
    /// DRAW command `lprops`: 1 uses in 1 converted scripts.
    ///
    /// Example: `lprops result`
    pub fn lprops(&mut self, args: &[&str]) -> DrawResult {
        // `lprops shape`: the length of its edges and their centre of mass.
        need(args, 1, "lprops shape")?;
        let shape = self.shape_var(args[0])?;
        Ok(mass_report(properties::length(&self.model, shape)))
    }
    /// DRAW command `mat`: 10 uses in 7 converted scripts.
    ///
    /// Example: `mat`
    pub fn mat(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("mat", args)
    }
    /// DRAW command `maxtolerance`: 13 uses in 8 converted scripts.
    ///
    /// Example: `maxtolerance a`
    pub fn maxtolerance(&mut self, args: &[&str]) -> DrawResult {
        // `maxtolerance shape`: the least and greatest tolerances of its faces, edges and
        // vertices, as OCCT reports them.
        need(args, 1, "maxtolerance shape")?;
        let shape = self.shape_var(args[0])?;
        let range = |kind: NodeKind| {
            let all = self.model.unique_subshapes(shape, kind);
            let tolerances: Vec<f64> = all.iter().map(|&s| self.model.tolerance(s)).collect();
            let min = tolerances.iter().copied().fold(f64::INFINITY, f64::min);
            let max = tolerances.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            (all.len(), min, max)
        };
        let (f, e, v) = (range(NodeKind::Face), range(NodeKind::Edge), range(NodeKind::Vertex));
        let mut out = format!("\n## Tolerances on the shape {}  (nbFaces:{}  nbEdges:{} nbVtx:{})\n", args[0], f.0, e.0, v.0);
        for (name, (_, min, max)) in [("Face  ", f), ("Edge  ", e), ("Vertex", v)] {
            if min <= max {
                out += &format!("\n    {name} : Min {}    Max  {} \n ", format_c_exponent(min), format_c_exponent(max));
            }
        }
        Ok(out)
    }
    /// DRAW command `middlepath`: 15 uses in 15 converted scripts.
    ///
    /// Example: `middlepath result b b_3 b_4`
    pub fn middlepath(&mut self, args: &[&str]) -> DrawResult {
        // `middlepath result shape start end`: the path through the middle of a pipe-like
        // shape from one of its sections to another.
        let [result, shape, start, end] = args else { return self.not_implemented("middlepath", args) };
        let (shape, start, end) = (self.shape_var(shape)?, self.shape_var(start)?, self.shape_var(end)?);
        let path = sweep::middle_path(&mut self.model, shape, start, end)?;
        self.set_shape_var(result, path)
    }
    /// DRAW command `mk2dcurve`: 7 uses in 5 converted scripts.
    ///
    /// Example: `mk2dcurve cc ac a`
    pub fn mk2dcurve(&mut self, args: &[&str]) -> DrawResult {
        // `mk2dcurve name edge [face | index]`: the edge's curve on the face (the index-th face
        // using it; the first if none given), over its range.
        need(args, 2, "mk2dcurve name edge [face|index]")?;
        let edge = self.shape_var(args[1])?;
        let edge = self.model.unique_subshapes(edge, NodeKind::Edge).into_iter().next().ok_or_else(|| DrawError::Failed(format!("{} has no edge", args[1])))?;
        let edge_id = kentron::topology::EdgeId::try_from(edge.node).expect("an edge");
        let face = match args.get(2) {
            Some(name) if matches!(self.get(name), Some(DrawObject::Shape(_))) => Some(self.shape_var(name)?),
            _ => None,
        };
        let coedge = match face {
            Some(f) => {
                let f = self.model.unique_subshapes(f, NodeKind::Face).into_iter().next().ok_or_else(|| DrawError::Failed("mk2dcurve: no face".into()))?;
                self.model.find_coedge(edge_id, kentron::topology::FaceId::try_from(f.node).expect("a face"))
            }
            None => {
                let index = match args.get(2) { Some(i) => self.number(i)? as usize, None => 1 };
                let coedges: Vec<kentron::topology::CoEdgeId> = self.model.parents(edge_id.into()).into_iter().filter_map(|n| kentron::topology::CoEdgeId::try_from(n).ok()).collect();
        coedges.get(index.saturating_sub(1)).copied()
            }
        };
        let coedge = coedge.ok_or_else(|| DrawError::Failed("mk2dcurve: the edge has no curve on that face".into()))?;
        let (c, first, last) = self.model.effective_pcurve(coedge).ok_or_else(|| DrawError::Failed("mk2dcurve: no pcurve".into()))?;
        self.set(args[0], DrawObject::Curve2(Curve2::Trimmed(TrimmedCurve2::new(c, first, last))));
        Ok(String::new())
    }
    /// DRAW command `mkcurve`: 35 uses in 24 converted scripts.
    ///
    /// Example: `mkcurve cu11 res1_1`
    pub fn mkcurve(&mut self, args: &[&str]) -> DrawResult {
        // `mkcurve name edge`: the edge's 3D curve over its range.
        need(args, 2, "mkcurve name edge")?;
        let edge = self.shape_var(args[1])?;
        let edge = self.model.unique_subshapes(edge, NodeKind::Edge).into_iter().next().ok_or_else(|| DrawError::Failed(format!("{} has no edge", args[1])))?;
        let (c, first, last) = self.model.edge_curve(Shape { orientation: Orientation::Forward, ..edge }).ok_or_else(|| DrawError::Failed("mkcurve: a degenerate edge".into()))?;
        self.set(args[0], DrawObject::Curve(Curve::Trimmed(TrimmedCurve::new(c, first, last))));
        Ok(String::new())
    }
    /// DRAW command `mkedge`: 432 uses in 257 converted scripts.
    ///
    /// Example: `mkedge w1 w1 0 pi`
    pub fn mkedge(&mut self, args: &[&str]) -> DrawResult {
        // `mkedge edge curve [first last]`: an edge on a 3D curve, over its whole domain or
        // the given parameter range.
        need(args, 2, "mkedge edge curve [first last]")?;
        let curve = self.curve_var(args[1])?;
        let edge = match args.len() {
            2 => make::edge_from_curve(&mut self.model, curve)?,
            4 => {
                let (first, last) = (self.number(args[2])?, self.number(args[3])?);
                make::edge_from_curve_range(&mut self.model, curve, first, last)?
            }
            _ => return self.not_implemented("mkedge", args),
        };
        self.set_shape_var(args[0], edge)
    }
    /// DRAW command `mkevol`: 103 uses in 103 converted scripts.
    ///
    /// Example: `mkevol result s`
    pub fn mkevol(&mut self, args: &[&str]) -> DrawResult {
        // `mkevol result shape`: starts an evolving fillet of the shape.
        need(args, 2, "mkevol result shape")?;
        let shape = self.shape_var(args[1])?;
        self.evolving = Some((args[0].to_string(), shape, Vec::new()));
        Ok(String::new())
    }
    /// DRAW command `mkface`: 785 uses in 305 converted scripts.
    ///
    /// Example: `mkface face surf wire`
    pub fn mkface(&mut self, args: &[&str]) -> DrawResult {
        // `mkface name surface [wire | umin umax vmin vmax]`: the face of the surface bounded
        // by the wire, or on the given rectangle of its parameters, or on its whole domain.
        need(args, 2, "mkface name surface [wire | umin umax vmin vmax]")?;
        let surface = self.surface_var(args[1])?;
        let face = match args.len() {
            2 => make::face_on_surface(&mut self.model, surface)?,
            3 => {
                let wire = self.shape_var(args[2])?;
                make::face_on(&mut self.model, surface, wire)?
            }
            6 => {
                let v = self.numbers(&args[2..])?;
                make::face_on_domain(&mut self.model, surface, v[0], v[1], v[2], v[3])?
            }
            _ => return self.not_implemented("mkface", args),
        };
        self.set_shape_var(args[0], face)
    }
    /// DRAW command `mkoffset`: 19 uses in 17 converted scripts.
    ///
    /// Example: `mkoffset resoffset a 1 20`
    pub fn mkoffset(&mut self, args: &[&str]) -> DrawResult {
        // `mkoffset result shape count step [a|i [alt]]`: offsets at step, 2 step, ... named
        // result_1, result_2, ...
        let [result, shape, count, step, rest @ ..] = args else { return self.not_implemented("mkoffset", args) };
        if !rest.is_empty() && rest[0] != "a" {
            return self.not_implemented("mkoffset", args);
        }
        let number = |a: &str| a.parse::<f64>().map_err(|_| DrawError::Failed(format!("mkoffset: {a} is not a number")));
        let (count, step) = (number(count)? as usize, number(step)?);
        let wire = self.shape_var(shape)?;
        for k in 1..=count {
            let offset = kentron::offset::planar_wire(&mut self.model, wire, step * k as f64, Default::default())?;
            self.set_shape_var(&format!("{result}_{k}"), offset)?;
        }
        Ok(String::new())
    }
    /// Not a DRAW command: the offset group's end script's counts of the result's vertices,
    /// edges and wires (an empty result to have no wires).
    pub fn offset_counts(&mut self, args: &[&str]) -> DrawResult {
        let [result, v, e, w] = args else { return self.not_implemented("offset_counts", args) };
        let w = if *v == "0" { "0" } else { w };
        self.checknbshapes(&[result, "-vertex", v, "-edge", e, "-wire", w])
    }
    /// Not a DRAW command: the offset group's end script, `mkoffset result s 1 [expr
    /// length(s) * fraction]`, with the result named `result`.
    pub fn mkoffset_fraction(&mut self, args: &[&str]) -> DrawResult {
        let [result, shape, fraction] = args else { return self.not_implemented("mkoffset_fraction", args) };
        let fraction = fraction.parse::<f64>().map_err(|_| DrawError::Failed(format!("mkoffset_fraction: {fraction}?")))?;
        let wire = self.shape_var(shape)?;
        let distance = kentron::properties::length(&self.model, wire).mass * fraction;
        let offset = kentron::offset::planar_wire(&mut self.model, wire, distance, Default::default())?;
        self.set_shape_var(result, offset)
    }
    /// DRAW command `mkoffsetshape`: 1 uses in 1 converted scripts.
    ///
    /// Example: `mkoffsetshape result sh 1 1e-7 1 0 i`
    pub fn mkoffsetshape(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("mkoffsetshape", args)
    }
    /// DRAW command `mkplane`: 509 uses in 362 converted scripts.
    ///
    /// Example: `mkplane w w`
    pub fn mkplane(&mut self, args: &[&str]) -> DrawResult {
        // `mkplane name wire [OnlyPlane]`: the planar face bounded by the wire.
        need(args, 2, "mkplane name wire [OnlyPlane]")?;
        let wire = self.shape_var(args[1])?;
        // (An open wire gives a face on its plane, as OCCT's face bounded by an open wire,
        // good for its surface.)
        let face = match make::planar_face(&mut self.model, wire) {
            Err(kentron::Error::InvalidInput(m)) if m.contains("closed wire") => {
                let plane = make::plane_of_wire(&self.model, wire)?;
                let b = properties::bounding_box(&self.model, wire);
                let size = (b.max - b.min).length();
                let (u, v) = plane.parameters(b.min + (b.max - b.min) * 0.5);
                make::face_on_domain(&mut self.model, Surface::Plane(plane), u - size, u + size, v - size, v + size)?
            }
            other => other?,
        };
        self.set_shape_var(args[0], face)
    }
    /// DRAW command `mkpoint`: 1 uses in 1 converted scripts.
    ///
    /// Example: `mkpoint p a_2`
    pub fn mkpoint(&mut self, args: &[&str]) -> DrawResult {
        // `mkpoint name vertex`: the vertex's point.
        need(args, 2, "mkpoint name vertex")?;
        let v = self.shape_var(args[1])?;
        let p = self.model.point(v);
        self.set(args[0], DrawObject::Point(p));
        Ok(String::new())
    }
    /// DRAW command `mksurface`: 89 uses in 75 converted scripts.
    ///
    /// Example: `mksurface surf a_26`
    pub fn mksurface(&mut self, args: &[&str]) -> DrawResult {
        // `mksurface name face`: the surface of a face (placed where the face is).
        need(args, 2, "mksurface name face")?;
        let face = self.shape_var(args[1])?;
        if face.kind() != NodeKind::Face {
            return Err(DrawError::Failed(format!("mksurface: {} is not a face", args[1])));
        }
        let surface = self.model.surface(face).cloned().ok_or_else(|| DrawError::Failed(format!("{} has no surface", args[1])))?;
        self.set(args[0], DrawObject::Surface(surface.transformed(&face.placement)));
        Ok(String::new())
    }
    /// DRAW command `mksweep`: 108 uses in 108 converted scripts.
    ///
    /// Example: `mksweep f_1`
    pub fn mksweep(&mut self, args: &[&str]) -> DrawResult {
        // `mksweep spine [-C [angle tolerance]]`: a sweep along the spine begins (OCCT's
        // BRepOffsetAPI_MakePipeShell); the corrected Frenet trihedron, no profiles yet.
        // (`-C` merges the spine's edges that carry on smoothly into one: kentron sweeps
        // along them as they are, smoothly from one to the next.)
        let args: &[&str] = match args {
            [spine, "-C", rest @ ..] if rest.len() <= 2 => std::slice::from_ref(spine),
            _ => args,
        };
        if args.len() != 1 {
            return self.not_implemented("mksweep", args);
        }
        let spine = self.shape_var(args[0])?;
        self.sweep_spine = Some(spine);
        self.sweep_mode = vec!["-CF".into()];
        self.sweep_profiles.clear();
        Ok(String::new())
    }
    /// DRAW command `mkvolume`: 99 uses in 99 converted scripts.
    ///
    /// Example: `mkvolume result f1 f2`
    pub fn mkvolume(&mut self, args: &[&str]) -> DrawResult {
        // `mkvolume result shapes... [-c] [-ni] [-ai]`: the solids the shapes' faces enclose
        // (`-c`: the compounds' parts intersected with each other; `-ni`: the faces as they
        // are, not intersected; `-ai`: no internal shapes kept, which kentron never keeps).
        let unsplit = args.contains(&"-ni");
        let expand = args.contains(&"-c");
        let names: Vec<&str> = args.iter().copied().filter(|a| !a.starts_with('-')).collect();
        need(&names, 2, "mkvolume result shapes... [-c] [-ni] [-ai]")?;
        let mut shapes = Vec::new();
        for name in &names[1..] {
            let shape = self.shape_var(name)?;
            if expand && shape.kind() == NodeKind::Compound {
                shapes.extend(self.model.children(shape));
            } else {
                shapes.push(shape);
            }
        }
        let result = if unsplit {
            kentron::boolean::enclosed_solids_unsplit(&mut self.model, &shapes)?
        } else {
            kentron::boolean::enclosed_solids(&mut self.model, &shapes)?
        };
        self.set_shape_var(names[0], result)
    }
    /// DRAW command `modified`: 9 uses in 5 converted scripts.
    ///
    /// Example: `modified f_mod h b2_12`
    pub fn modified(&mut self, args: &[&str]) -> DrawResult {
        // `modified result history shape`: the shapes `shape` was modified into (one alone,
        // else a compound of them).
        need(args, 3, "modified result history shape")?;
        let history = self.draw_history(args[1])?;
        let input = self.shape_var(args[2])?;
        let shapes: Vec<Shape> = history.modified.iter().filter(|(node, _)| *node == input.node).flat_map(|(_, s)| s.clone()).collect();
        match shapes[..] {
            [] => Ok("The shape has not been modified\n".into()),
            [one] => self.set_shape_var(args[0], one),
            _ => {
                let compound = self.compound_of(&shapes)?;
                self.set_shape_var(args[0], compound)
            }
        }
    }
    /// DRAW command `isdeleted`.
    ///
    /// Example: `isdeleted rf_hist s_2`
    pub fn isdeleted(&mut self, args: &[&str]) -> DrawResult {
        // `isdeleted history shape`: whether the operation deleted `shape`.
        need(args, 2, "isdeleted history shape")?;
        let history = self.draw_history(args[0])?;
        let input = self.shape_var(args[1])?;
        Ok(if history.deleted.contains(&input.node) { "Deleted." } else { "Not deleted." }.into())
    }
    /// DRAW command `removefeatures`.
    ///
    /// Example: `removefeatures result s s_2 s_3`
    pub fn removefeatures(&mut self, args: &[&str]) -> DrawResult {
        // `removefeatures result shape features... [-parallel]`: the shape without the faces
        // of the features (faces, or compounds of them).
        let names: Vec<&str> = args.iter().copied().filter(|a| !a.starts_with('-')).collect();
        if names.len() < 3 {
            return self.not_implemented("removefeatures", args);
        }
        let shape = self.shape_var(names[1])?;
        let mut faces = Vec::new();
        for name in &names[2..] {
            let feature = self.shape_var(name)?;
            faces.extend(self.model.unique_subshapes(feature, NodeKind::Face));
        }
        let removed = kentron::defeature::remove_faces(&mut self.model, shape, &faces)?;
        let shapes_of = |model: &kentron::topology::Model, nodes: &[kentron::topology::NodeId]| -> Vec<Shape> {
            nodes.iter().filter_map(|&n| model.unique_subshapes(removed.shape, n.kind).into_iter().find(|s| s.node == n)).collect()
        };
        let mut modified = Vec::new();
        let mut deleted = Vec::new();
        for kind in [NodeKind::Face, NodeKind::Edge, NodeKind::Vertex] {
            for s in self.model.unique_subshapes(shape, kind) {
                if removed.history.is_deleted(s.node) {
                    deleted.push(s.node);
                } else if let Some(images) = removed.history.modified(s.node) {
                    modified.push((s.node, shapes_of(&self.model, images)));
                }
            }
        }
        self.last_history = DrawHistory { modified, deleted, ..Default::default() };
        self.set_shape_var(names[0], removed.shape)
    }
    /// DRAW command `movep`: 4 uses in 4 converted scripts.
    ///
    /// Example: `movep p 2 2 0 0 1`
    pub fn movep(&mut self, args: &[&str]) -> DrawResult {
        // `movep surface row col dx dy dz`: the B-spline or Bézier surface's pole (from 1) moved.
        need(args, 6, "movep surface row col dx dy dz")?;
        let (i, j) = (self.number(args[1])? as usize, self.number(args[2])? as usize);
        let d = self.numbers(&args[3..6])?;
        let d = Vector::new(d[0], d[1], d[2]);
        let s = match self.surface_var(args[0])? {
            Surface::BSpline(mut b) if i >= 1 && j >= 1 => { let p = b.pole(i - 1, j - 1); b.set_pole(i - 1, j - 1, p + d); Surface::BSpline(b) }
            Surface::Bezier(mut b) if i >= 1 && j >= 1 => { let p = b.pole(i - 1, j - 1); b.set_pole(i - 1, j - 1, p + d); Surface::Bezier(b) }
            _ => return Err(DrawError::Failed(format!("movep: {} is not a B-spline or Bézier surface", args[0]))),
        };
        self.set(args[0], DrawObject::Surface(s));
        Ok(String::new())
    }
    /// DRAW command `mu4`: 1 uses in 1 converted scripts.
    ///
    /// Example: `mu4`
    pub fn mu4(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("mu4", args)
    }
    /// DRAW command `nbshapes`: 18 uses in 14 converted scripts.
    ///
    /// Example: `nbshapes b`
    pub fn nbshapes(&mut self, args: &[&str]) -> DrawResult {
        // `nbshapes shape [-t]`: the numbers of its distinct sub-shapes of each kind, in
        // OCCT's report: each shared shape once however placed (each placement once with
        // `-t`).
        need(args, 1, "nbshapes shape [-t]")?;
        let shape = self.shape_var(args[0])?;
        let shared = !args.contains(&"-t");
        let mut out = format!("Number of shapes in {}\n", args[0]);
        let mut total = 0;
        for (name, kind) in SHAPE_KINDS {
            let n = count_subshapes(&self.model, shape, kind, shared);
            total += n;
            out += &format!(" {:<10}: {n}\n", name.to_uppercase());
        }
        out += &format!(" {:<10}: {total}\n", "SHAPE");
        Ok(out)
    }
    /// DRAW command `nexplode`: 328 uses in 306 converted scripts.
    ///
    /// Example: `nexplode s e`
    pub fn nexplode(&mut self, args: &[&str]) -> DrawResult {
        // `nexplode name F|E|V`: like `explode`, ordered by position: by 999 x + 99 y + 0.9 z
        // of each one's centre (of its edges' length for faces and edges), stably.
        need(args, 2, "nexplode name F|E|V")?;
        let shape = self.shape_var(args[0])?;
        let kind = match args[1].to_ascii_uppercase().as_str() {
            "F" | "FACE" => NodeKind::Face,
            "E" | "EDGE" => NodeKind::Edge,
            "V" | "VERTEX" => NodeKind::Vertex,
            t => return Err(DrawError::Failed(format!("nexplode: unknown type {t}"))),
        };
        let mut parts: Vec<Shape> = Vec::new();
        for s in self.model.explore(shape, kind) {
            if !parts.iter().any(|q| q.is_same(&s)) {
                parts.push(s);
            }
        }
        let key = |model: &kentron::topology::Model, s: Shape| {
            let p = if kind == NodeKind::Vertex { model.point(s) } else { properties::length(model, s).center_of_mass };
            p.x * 999.0 + p.y * 99.0 + p.z * 0.9
        };
        let mut keyed: Vec<(f64, Shape)> = parts.into_iter().map(|s| (key(&self.model, s), s)).collect();
        keyed.sort_by(|a, b| a.0.total_cmp(&b.0));
        let names: Vec<String> = (1..=keyed.len()).map(|i| format!("{}_{i}", args[0])).collect();
        for (name, (_, part)) in names.iter().zip(keyed) {
            self.set_shape_var(name, part)?;
        }
        Ok(names.join(" "))
    }
    /// DRAW command `normals`: 2 uses in 2 converted scripts.
    ///
    /// Example: `normals result`
    pub fn normals(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("normals", args)
    }
    /// DRAW command `nproject`: 21 uses in 21 converted scripts.
    ///
    /// Example: `nproject result w shell`
    pub fn nproject(&mut self, args: &[&str]) -> DrawResult {
        // `nproject result shape ... target [-g] [-d [dmax]] [tol [continuity [maxdeg [maxseg]]]]`:
        // the edges of the shapes projected along the normals of the target's faces, within
        // the faces (anywhere on their surfaces with -g) and dmax (1e-3) of them with -d.
        // The fits' continuity, degree and segments are kentron's own.
        need(args, 3, "nproject result shape ... target [-g] [-d [dmax]] [tol ...]")?;
        let mut shapes = Vec::new();
        let mut i = 1;
        while let Some(s) = args.get(i).and_then(|a| self.shape_var(a).ok()) {
            shapes.push(s);
            i += 1;
        }
        let Some(target) = shapes.pop().filter(|_| !shapes.is_empty()) else {
            return Err(DrawError::Failed("nproject: needs shapes and a target".into()));
        };
        let mut options = kentron::project::NormalProjection::default();
        if args.get(i) == Some(&"-g") {
            options.within_faces = false;
            i += 1;
        }
        if args.get(i) == Some(&"-d") {
            i += 1;
            let mut distance = 1e-3;
            if let Some(d) = args.get(i) {
                distance = self.number(d)?;
                i += 1;
            }
            options.max_distance = Some(distance);
        }
        if let Some(tol) = args.get(i) {
            options.tolerance = self.number(tol)?.max(1e-10);
        }
        let r = kentron::project::normal(&mut self.model, &shapes, target, &options)?;
        self.set_shape_var(args[0], r)
    }
    /// DRAW command `numshapes`: 2 uses in 2 converted scripts.
    ///
    /// Example: `numshapes shell`
    pub fn numshapes(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("numshapes", args)
    }
    /// DRAW command `nurbsconvert`: 119 uses in 113 converted scripts.
    ///
    /// Example: `nurbsconvert b1 b1`
    pub fn nurbsconvert(&mut self, args: &[&str]) -> DrawResult {
        // `nurbsconvert result shape [result shape ...]`: copies with B-spline geometry.
        if args.is_empty() || args.len() % 2 != 0 {
            return Err(DrawError::Failed("usage: nurbsconvert result shape [result shape ...]".into()));
        }
        for pair in args.chunks(2) {
            let shape = self.shape_var(pair[1])?;
            let converted = convert::to_nurbs(&mut self.model, shape)?;
            self.set_shape_var(pair[0], converted)?;
        }
        Ok(String::new())
    }
    /// DRAW command `offset`: 128 uses in 87 converted scripts.
    ///
    /// Example: `offset o1 ss1 -0.1`
    pub fn offset(&mut self, args: &[&str]) -> DrawResult {
        // `offset name basis distance [dx dy dz]`: the offset of a 2D curve, of a 3D curve
        // (in the plane normal to the direction) or of a surface.
        need(args, 3, "offset name basis distance [dx dy dz]")?;
        let distance = self.number(args[2])?;
        let object = match self.get(args[1]).cloned() {
            Some(DrawObject::Curve2(c)) => DrawObject::Curve2(Curve2::Offset(kentron::geometry::curve2::OffsetCurve2::new(c, distance))),
            Some(DrawObject::Curve(c)) if args.len() == 6 => {
                let v = self.numbers(&args[3..6])?;
                let d = Direction::try_from_vector(Vector::new(v[0], v[1], v[2])).ok_or_else(|| DrawError::Failed("offset: null direction".into()))?;
                DrawObject::Curve(Curve::Offset(kentron::geometry::curve::OffsetCurve::new(c, distance, d)))
            }
            Some(DrawObject::Surface(s)) => DrawObject::Surface(Surface::Offset(kentron::geometry::surface::OffsetSurface::new(s, distance))),
            _ => return self.not_implemented("offset", args),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `offset2dcurve`: 2 uses in 2 converted scripts.
    ///
    /// Example: `offset2dcurve o cc .5`
    pub fn offset2dcurve(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("offset2dcurve", args)
    }
    /// DRAW command `offsetcompshape`: 7 uses in 7 converted scripts.
    ///
    /// Example: `offsetcompshape d_thick draft -5 draft_2`
    pub fn offsetcompshape(&mut self, args: &[&str]) -> DrawResult {
        // `offsetcompshape result shape offset [tol] [face ...]`: as `offsetshape`, the faces'
        // offsets intersected completely (as kentron's always are).
        self.offset_shape_command(args, kentron::offset::Join::Arc)
    }
    /// DRAW command `offsetcurve`: 2 uses in 2 converted scripts.
    ///
    /// Example: `offsetcurve o cc .5 pp`
    pub fn offsetcurve(&mut self, args: &[&str]) -> DrawResult {
        // `offsetcurve result curve distance direction`: the curve offset by the distance,
        // square to it and to the direction (a point variable read as a vector), as OCCT's
        // Geom_OffsetCurve.
        let [result, curve, distance, direction] = args else { return self.not_implemented("offsetcurve", args) };
        let basis = self.curve_var(curve)?;
        let distance = self.number(distance)?;
        let Some(DrawObject::Point(p)) = self.get(direction).cloned() else {
            return Err(DrawError::Failed(format!("offsetcurve: {direction} is not a point")));
        };
        let direction = kentron::math::Direction::try_from_vector(kentron::math::Vector::new(p.x, p.y, p.z))
            .ok_or_else(|| DrawError::Failed("offsetcurve: a zero direction".into()))?;
        let offset = kentron::geometry::curve::Curve::Offset(kentron::geometry::curve::OffsetCurve::new(basis, distance, direction));
        self.set(result, DrawObject::Curve(offset));
        Ok(String::new())
    }
    /// DRAW command `offsetload`: 174 uses in 174 converted scripts.
    ///
    /// Example: `offsetload s 0.14`
    pub fn offsetload(&mut self, args: &[&str]) -> DrawResult {
        // `offsetload shape distance [faces...]`: faces left out make a thick solid.
        let [shape, distance, faces @ ..] = args else { return self.not_implemented("offsetload", args) };
        let distance = distance.parse::<f64>().map_err(|_| DrawError::Failed(format!("offsetload: {distance}?")))?;
        let shape = self.shape_var(shape)?;
        let faces = faces.iter().map(|f| self.shape_var(f)).collect::<Result<Vec<_>, _>>()?;
        self.offset_load = Some((shape, distance, faces));
        self.offset_on_faces.clear();
        Ok(String::new())
    }
    /// DRAW command `offsetonface`: 290 uses in 49 converted scripts.
    ///
    /// Example: `offsetonface b1_1 1.0`
    pub fn offsetonface(&mut self, args: &[&str]) -> DrawResult {
        // `offsetonface face distance [face distance ...]`: faces of the loaded shape moved by
        // distances of their own.
        if args.is_empty() || args.len() % 2 != 0 {
            return self.not_implemented("offsetonface", args);
        }
        for pair in args.chunks(2) {
            let face = self.shape_var(pair[0])?;
            let d = self.number(pair[1])?;
            self.offset_on_faces.push((face, d));
        }
        Ok(String::new())
    }
    /// DRAW command `offsetparameter`: 174 uses in 174 converted scripts.
    ///
    /// Example: `offsetparameter 1.e-7 c i`
    pub fn offsetparameter(&mut self, args: &[&str]) -> DrawResult {
        // `offsetparameter tol c|p a|i [r|u]`: complete or partial intersection (kentron's
        // offsets find all intersections either way), arc or intersection joins; with `r`,
        // faces split by edges inside the result's planes and surfaces are merged (OCCT's
        // RemoveIntEdges).
        match args {
            [_, "c" | "p", join, rest @ ..] if rest.is_empty() || rest == ["u"] || rest == ["r"] => {
                self.offset_merge = rest == ["r"];
                self.offset_join = match *join {
                    "a" => kentron::offset::Join::Arc,
                    "i" => kentron::offset::Join::Intersection,
                    _ => return self.not_implemented("offsetparameter", args),
                };
                Ok(String::new())
            }
            _ => self.not_implemented("offsetparameter", args),
        }
    }
    /// DRAW command `offsetperform`: 174 uses in 174 converted scripts.
    ///
    /// Example: `offsetperform result`
    pub fn offsetperform(&mut self, args: &[&str]) -> DrawResult {
        let [result] = args else { return self.not_implemented("offsetperform", args) };
        let Some((shape, distance, faces)) = self.offset_load.clone() else {
            return Err(DrawError::Failed("offsetperform: no offsetload".into()));
        };
        let join = self.offset_join;
        let r = if faces.is_empty() {
            kentron::offset::shape_with_faces(&mut self.model, shape, distance, &self.offset_on_faces.clone(), join)?
        } else if !self.offset_on_faces.is_empty() {
            return self.not_implemented("offsetperform", args);
        } else {
            kentron::offset::thick_solid(&mut self.model, shape, &faces, distance, join)?
        };
        let r = if self.offset_merge { kentron::unify::unify_same_domain(&mut self.model, r)? } else { r };
        self.set_shape_var(result, r)
    }
    /// DRAW command `offsetshape`: 9 uses in 9 converted scripts.
    ///
    /// Example: `offsetshape result a 10`
    pub fn offsetshape(&mut self, args: &[&str]) -> DrawResult {
        // `offsetshape result shape offset [tol] [face ...]`: the shape offset with arcs
        // joining its faces; a thick solid with the faces given left open.
        self.offset_shape_command(args, kentron::offset::Join::Arc)
    }

    /// `offsetshape`, `offsetcompshape` and `offsetinter` (kentron's offsets always
    /// intersect the faces completely), with the given join.
    fn offset_shape_command(&mut self, args: &[&str], join: kentron::offset::Join) -> DrawResult {
        need(args, 3, "offsetshape result shape offset [tol] [face ...]")?;
        let shape = self.shape_var(args[1])?;
        let distance: f64 = args[2].parse().map_err(|_| DrawError::Failed(format!("offsetshape: bad offset {}", args[2])))?;
        let mut faces = Vec::new();
        for (k, name) in args[3..].iter().enumerate() {
            match self.shape_var(name) {
                Ok(face) if face.kind() == NodeKind::Face => faces.push(face),
                // (A tolerance first.)
                _ if k == 0 && name.parse::<f64>().is_ok() => {}
                _ => return Err(DrawError::Failed(format!("offsetshape: {name} is not a face"))),
            }
        }
        let r = if faces.is_empty() {
            kentron::offset::shape(&mut self.model, shape, distance, join)?
        } else {
            kentron::offset::thick_solid(&mut self.model, shape, &faces, distance, join)?
        };
        self.set_shape_var(args[0], r)
    }
    /// DRAW command `offsetshapesimple`: 2 uses in 1 converted scripts.
    ///
    /// Example: `offsetshapesimple res1 fix1 10.`
    pub fn offsetshapesimple(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("offsetshapesimple", args)
    }
    /// DRAW command `openoffset`: 32 uses in 31 converted scripts.
    ///
    /// Example: `openoffset res1 ww 1 12.5`
    pub fn openoffset(&mut self, args: &[&str]) -> DrawResult {
        // `openoffset result wire count step [a|i] [-approx]`: the open wire's offsets on one
        // side at step, 2 step, ..., named result_1, result_2, ... (`-approx`, approximating
        // the wire by lines and arcs first, is not needed: its curves are offset as they are).
        let args: Vec<&str> = args.iter().copied().filter(|a| *a != "-approx").collect();
        let [result, shape, count, step, rest @ ..] = &args[..] else { return self.not_implemented("openoffset", &args) };
        let join = match rest {
            [] | ["a"] => kentron::offset::Join::Arc,
            ["i"] => kentron::offset::Join::Intersection,
            _ => return self.not_implemented("openoffset", &args),
        };
        let (count, step) = (self.number(count)? as usize, self.number(step)?);
        let wire = self.shape_var(shape)?;
        for k in 1..=count {
            let offset = kentron::offset::open_planar_wire(&mut self.model, wire, step * k as f64, join)?;
            self.set_shape_var(&format!("{result}_{k}"), offset)?;
        }
        Ok(String::new())
    }
    /// DRAW command `orientation`: 49 uses in 31 converted scripts.
    ///
    /// Example: `orientation w2 R`
    pub fn orientation(&mut self, args: &[&str]) -> DrawResult {
        // `orientation name ... F|R|E|I`: sets the shapes' orientation.
        need(args, 2, "orientation name ... F|R|E|I")?;
        let n = args.len();
        let orientation = match args[n - 1].to_ascii_uppercase().as_str() {
            "F" => Orientation::Forward,
            "R" => Orientation::Reversed,
            "E" => Orientation::External,
            "I" => Orientation::Internal,
            o => return Err(DrawError::Failed(format!("orientation: unknown orientation {o}"))),
        };
        for name in &args[..n - 1] {
            let shape = self.shape_var(name)?;
            self.set_shape_var(name, Shape { orientation, ..shape })?;
        }
        Ok(String::new())
    }
    /// DRAW command `parabola`: 52 uses in 51 converted scripts.
    ///
    /// Example: `parabola w1 0 0 0 1 0 0 0 0 1 50`
    pub fn parabola(&mut self, args: &[&str]) -> DrawResult {
        // `parabola name x y [dx dy] focal` (2D) or `parabola name x y z [nx ny nz [ux uy uz]]
        // focal`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, f) = v.split_at(v.len().saturating_sub(1));
        let object = match frame_args.len() {
            2 | 4 => DrawObject::Curve2(Curve2::Parabola(Parabola2::new(Draw::curve_frame2(frame_args)?, f[0]))),
            _ => DrawObject::Curve(Curve::Parabola(Parabola::new(Draw::curve_frame(frame_args)?, f[0]))),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `parameters`: 2 uses in 1 converted scripts.
    ///
    /// Example: `parameters c x y z 1.e-7 u`
    pub fn parameters(&mut self, args: &[&str]) -> DrawResult {
        // `parameters surface x y z tolerance u v`, `parameters curve x y z tolerance u` or
        // `parameters curve2d x y tolerance u`: the parameters of a point on it, into
        // variables; an error if the point is farther than the tolerance.
        use kentron::geometry::extrema;
        let wrong = || DrawError::Failed("Wrong point".into());
        match (args.len(), self.get(args.first().copied().unwrap_or_default()).cloned()) {
            (7, Some(DrawObject::Surface(s))) => {
                let v = self.numbers(&args[1..5])?;
                let (u, w, d) = extrema::nearest_on_surface(Point::new(v[0], v[1], v[2]), &s, None);
                self.set(args[5], DrawObject::Number(u));
                self.set(args[6], DrawObject::Number(w));
                if d > v[3] { Err(wrong()) } else { Ok(String::new()) }
            }
            (6, Some(DrawObject::Curve(c))) => {
                let v = self.numbers(&args[1..5])?;
                let found = extrema::nearest_on_curve(Point::new(v[0], v[1], v[2]), &c, c.first_parameter(), c.last_parameter());
                self.set(args[5], DrawObject::Number(found.parameter));
                if found.distance > v[3] { Err(wrong()) } else { Ok(String::new()) }
            }
            (5, Some(DrawObject::Curve2(c))) => {
                let v = self.numbers(&args[1..4])?;
                let found = extrema::point_curve2(Point2::new(v[0], v[1]), &c).into_iter().min_by(|a, b| a.distance.total_cmp(&b.distance)).ok_or_else(wrong)?;
                self.set(args[4], DrawObject::Number(found.parameter));
                if found.distance > v[2] { Err(wrong()) } else { Ok(String::new()) }
            }
            _ => Err(DrawError::Failed("usage: parameters surface|curve x y [z] tolerance u [v]".into())),
        }
    }
    /// DRAW command `pbsplinecurve`: 13 uses in 13 converted scripts.
    ///
    /// Example: `pbsplinecurve sp1 3 9 0 2 73.198335334976 1 109.22594821708 1 168.29694729401 1 244.58155163942 1 30`
    pub fn pbsplinecurve(&mut self, args: &[&str]) -> DrawResult {
        // `pbsplinecurve name degree nknots k1 m1 ... x1 y1 z1 w1 ...`: a periodic B-spline.
        let (degree, knots, mults, rest) = self.bspline_knots(args, "pbsplinecurve")?;
        if rest.len() % 4 != 0 {
            return Err(DrawError::Failed("pbsplinecurve: poles come as x y z w".into()));
        }
        let poles: Vec<Point> = rest.chunks(4).map(|c| Point::new(c[0], c[1], c[2])).collect();
        let weights: Vec<f64> = rest.chunks(4).map(|c| c[3]).collect();
        let mut curve = BSplineCurve::try_new_periodic(poles, knots, mults, degree)?;
        if weights.iter().any(|&w| w != 1.0) {
            curve = curve.with_weights(weights);
        }
        self.set(args[0], DrawObject::Curve(Curve::BSpline(curve)));
        Ok(String::new())
    }
    /// DRAW command `pcone`: 386 uses in 295 converted scripts.
    ///
    /// Example: `pcone s 90 40 150 90`
    pub fn pcone(&mut self, args: &[&str]) -> DrawResult {
        // `pcone name [plane] r1 r2 height [angle]` (angle in degrees).
        need(args, 4, "pcone name [plane] r1 r2 height [angle]")?;
        let (frame, values) = self.primitive_frame(&args[1..])?;
        let cone = match values[..] {
            [r1, r2, h] => Cone::new(r1, r2, h),
            [r1, r2, h, angle] => Cone::new(r1, r2, h).with_angle(angle.to_radians()),
            _ => return Err(DrawError::Failed("pcone: syntax error".into())),
        };
        let shape = cone.with_frame(frame).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `pcurve`: 19 uses in 17 converted scripts.
    ///
    /// Example: `pcurve ff`
    pub fn pcurve(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("pcurve", args)
    }
    /// DRAW command `pcylinder`: 1140 uses in 931 converted scripts.
    ///
    /// Example: `pcylinder s 50 100 90`
    pub fn pcylinder(&mut self, args: &[&str]) -> DrawResult {
        // `pcylinder name [plane] radius height [angle]` (angle in degrees).
        need(args, 3, "pcylinder name [plane] radius height [angle]")?;
        let (frame, values) = self.primitive_frame(&args[1..])?;
        let cylinder = match values[..] {
            [r, h] => Cylinder::new(r, h),
            [r, h, angle] => Cylinder::new(r, h).with_angle(angle.to_radians()),
            _ => return Err(DrawError::Failed("pcylinder: syntax error".into())),
        };
        let shape = cylinder.with_frame(frame).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `pd`: 6 uses in 2 converted scripts.
    ///
    /// Example: `pd`
    pub fn pd(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("pd", args)
    }
    /// DRAW command `pipe`: 135 uses in 123 converted scripts.
    ///
    /// Example: `pipe s w profile`
    pub fn pipe(&mut self, args: &[&str]) -> DrawResult {
        // `pipe name spine profile [mode]`: the profile swept along the spine (mode 0, the
        // default corrected Frenet trihedron, is how kentron sweeps).
        need(args, 3, "pipe name spine profile")?;
        if args.len() > 4 || (args.len() == 4 && args[3] != "0") {
            return self.not_implemented("pipe", args);
        }
        let (spine, profile) = (self.shape_var(args[1])?, self.shape_var(args[2])?);
        let shape = sweep::pipe(&mut self.model, profile, spine)?.shape;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `pl`: 28 uses in 10 converted scripts.
    ///
    /// Example: `pl`
    pub fn pl(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("pl", args)
    }
    /// DRAW command `plane`: 698 uses in 265 converted scripts.
    ///
    /// Example: `plane pl1 100 100 100 0 0 1 1 0 0`
    pub fn plane(&mut self, args: &[&str]) -> DrawResult {
        // `plane name [x y z [nx ny nz [dx dy dz]]]`: origin, normal and X direction.
        need(args, 1, "plane name [x y z [nx ny nz [dx dy dz]]]")?;
        let values = args[1..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let vector = |i: usize| Vector::new(values[i], values[i + 1], values[i + 2]);
        let direction = |i: usize| {
            Direction::try_from_vector(vector(i)).ok_or_else(|| DrawError::Failed("plane: null direction".into()))
        };
        let frame = match values.len() {
            0 => Frame::WORLD,
            3 => Frame::from_z(Point::ORIGIN + vector(0), Direction::Z),
            6 => Frame::from_z(Point::ORIGIN + vector(0), direction(3)?),
            9 => Frame::new(Point::ORIGIN + vector(0), direction(3)?, direction(6)?),
            _ => return self.not_implemented("plane", args),
        };
        self.set(args[0], DrawObject::Surface(Surface::Plane(Plane::from_frame(frame))));
        Ok(String::new())
    }
    /// DRAW command `pmirror`: 6 uses in 6 converted scripts.
    ///
    /// Example: `pmirror result 0 0 0`
    pub fn pmirror(&mut self, args: &[&str]) -> DrawResult {
        // `pmirror name ... x y z`: the point symmetry.
        need(args, 4, "pmirror name ... x y z")?;
        let n = args.len();
        let v = self.numbers(&args[n - 3..])?;
        self.transform_vars(&args[..n - 3], &Transform::mirror_point(Point::new(v[0], v[1], v[2])))
    }
    /// DRAW command `point`: 52 uses in 15 converted scripts.
    ///
    /// Example: `point pnt1 0 -0.3 0`
    pub fn point(&mut self, args: &[&str]) -> DrawResult {
        // `point name x y [z]`: a 2D or 3D point.
        need(args, 3, "point name x y [z]")?;
        let v = self.numbers(&args[1..])?;
        let object = match v[..] {
            [x, y] => DrawObject::Point2(Point2::new(x, y)),
            [x, y, z] => DrawObject::Point(Point::new(x, y, z)),
            _ => return Err(DrawError::Failed("point: syntax error".into())),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `polygon2d`: 1 uses in 1 converted scripts.
    ///
    /// Example: `polygon2d result 3 0 0 1 0 0 1`
    pub fn polygon2d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("polygon2d", args)
    }
    /// DRAW command `polygonprops`: 2 uses in 1 converted scripts.
    ///
    /// Example: `polygonprops result -area anArea1 -perimeter aPerimeter1`
    pub fn polygonprops(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("polygonprops", args)
    }
    /// DRAW command `polyline`: 246 uses in 195 converted scripts.
    ///
    /// Example: `polyline w2 -15 0 0 0 -5 0 15 0 0`
    pub fn polyline(&mut self, args: &[&str]) -> DrawResult {
        // `polyline name x1 y1 z1 x2 y2 z2 ...`: straight edges through the points; closed if
        // the last point is the first.
        need(args, 7, "polyline name x1 y1 z1 x2 y2 z2 ...")?;
        if (args.len() - 1) % 3 != 0 {
            return Err(DrawError::Failed("polyline: coordinates come in threes".into()));
        }
        let points = args[1..]
            .chunks(3)
            .map(|c| Ok(Point::new(self.number(c[0])?, self.number(c[1])?, self.number(c[2])?)))
            .collect::<DrawResult<Vec<Point>>>()?;
        let closed = points.len() > 2 && points[0].distance(*points.last().unwrap()) <= tolerance::LINEAR;
        let distinct = if closed { &points[..points.len() - 1] } else { &points[..] };
        let vertices: Vec<_> = distinct.iter().map(|&p| make::vertex(&mut self.model, p)).collect();
        let count = if closed { vertices.len() } else { vertices.len() - 1 };
        let edges = (0..count)
            .map(|i| make::edge(&mut self.model, vertices[i], vertices[(i + 1) % vertices.len()]))
            .collect::<kentron::Result<Vec<_>>>()?;
        let wire = make::wire(&mut self.model, &edges)?;
        self.set_shape_var(args[0], wire)
    }
    /// DRAW command `polyvertex`: 22 uses in 19 converted scripts.
    ///
    /// Example: `polyvertex w v1 v2 v3 v1`
    pub fn polyvertex(&mut self, args: &[&str]) -> DrawResult {
        // `polyvertex name v1 v2 ...`: a wire of straight edges joining the vertices.
        need(args, 3, "polyvertex name v1 v2 ...")?;
        let vertices = args[1..].iter().map(|a| self.shape_var(a)).collect::<DrawResult<Vec<Shape>>>()?;
        let mut edges = Vec::new();
        for pair in vertices.windows(2) {
            edges.push(make::edge(&mut self.model, pair[0], pair[1])?);
        }
        let wire = make::wire(&mut self.model, &edges)?;
        self.set_shape_var(args[0], wire)
    }
    /// DRAW command `prism`: 485 uses in 317 converted scripts.
    ///
    /// Example: `prism s w 2 0 10`
    pub fn prism(&mut self, args: &[&str]) -> DrawResult {
        // `prism result base dx dy dz [Copy | Inf | SemiInf]`: the copy flag changes nothing in
        // kentron; infinite prisms are not supported.
        need(args, 5, "prism result base dx dy dz [Copy | Inf | SemiInf]")?;
        if args.len() > 5 && !args[5].eq_ignore_ascii_case("copy") {
            return self.not_implemented("prism", args);
        }
        let base = self.shape_var(args[1])?;
        let v = Vector::new(self.number(args[2])?, self.number(args[3])?, self.number(args[4])?);
        let result = sweep::extrude(&mut self.model, base, v)?;
        self.set_shape_var(args[0], result.shape)
    }
    /// DRAW command `prj`: 81 uses in 81 converted scripts.
    ///
    /// Example: `prj res wire plate 0 0 -1`
    pub fn prj(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("prj", args)
    }
    /// DRAW command `profile`: 406 uses in 290 converted scripts.
    ///
    /// Example: `profile rev S b_4 F 50 20 Y 50 C 10 180 Y -50 C 10 180`
    pub fn profile(&mut self, args: &[&str]) -> DrawResult {
        // `profile name [code values] ...`: DRAW's codes for the steps of `sketch::profile`.
        need(args, 1, "profile name [code values] ...")?;
        let mut steps = Vec::new();
        let mut outcome = Outcome::Face;
        let mut i = 1;
        while i < args.len() {
            // DRAW reads a code's first letter, and the second for XX, YY, TT, RR, IX, IY, WW.
            let letters = args[i].to_ascii_uppercase().into_bytes();
            let code = match (letters[0], letters.get(1).copied()) {
                (c @ (b'X' | b'Y' | b'T' | b'R' | b'W'), Some(d)) if c == d => format!("{}{}", c as char, c as char),
                (b'I', Some(d @ (b'X' | b'Y'))) => format!("I{}", d as char),
                (c, _) => (c as char).to_string(),
            };
            let count = match code.as_str() {
                "F" | "T" | "TT" | "D" | "C" => 2,
                "O" => 3,
                "P" => 6,
                "X" | "XX" | "Y" | "YY" | "L" | "R" | "RR" | "IX" | "IY" | "S" => 1,
                "W" | "WW" => 0,
                _ => return Err(DrawError::Failed(format!("profile : unknown code {}", args[i]))),
            };
            let values = &args[i + 1..(i + 1 + count).min(args.len())];
            if values.len() < count {
                return Err(DrawError::Failed("profile : bad number of arguments".into()));
            }
            i += 1 + count;
            if code == "S" {
                let shape = self.shape_var(values[0])?;
                match self.model.surface(shape) {
                    Some(Surface::Plane(p)) if shape.kind() == NodeKind::Face => {
                        steps.push(Step::Frame(*p.transformed(&shape.placement).frame()));
                    }
                    _ => return self.not_implemented("profile", args),
                }
                continue;
            }
            let v = values.iter().map(|x| self.number(x)).collect::<DrawResult<Vec<f64>>>()?;
            steps.push(match code.as_str() {
                "F" => Step::Start(v[0], v[1]),
                "O" => Step::Origin(Point::new(v[0], v[1], v[2])),
                "P" => Step::Plane { normal: Vector::new(v[0], v[1], v[2]), x: Vector::new(v[3], v[4], v[5]) },
                "X" => Step::AlongX(v[0]),
                "XX" => Step::ToX(v[0]),
                "Y" => Step::AlongY(v[0]),
                "YY" => Step::ToY(v[0]),
                "L" => Step::Forward(v[0]),
                "T" => Step::By(v[0], v[1]),
                "TT" => Step::To(v[0], v[1]),
                "R" => Step::Turn(v[0]),
                "RR" => Step::Heading(v[0]),
                "D" => Step::Direction(v[0], v[1]),
                "C" => Step::Arc { radius: v[0], angle: v[1] },
                "IX" => Step::ForwardToX(v[0]),
                "IY" => Step::ForwardToY(v[0]),
                // The rest of the codes are ignored.
                _ => {
                    outcome = if code == "W" { Outcome::ClosedWire } else { Outcome::OpenWire };
                    break;
                }
            });
        }
        let shape = sketch::profile(&mut self.model, Frame::WORLD, &steps, outcome)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `proj`: 5 uses in 5 converted scripts.
    ///
    /// Example: `proj s 12.0565, 36.7424, -9.68517`
    pub fn proj(&mut self, args: &[&str]) -> DrawResult {
        // `proj curve|surface x y z`: the point's projections onto the curve (all of them) or the
        // surface (the nearest): each as `ext_i`, the segment from the point to it (the projected
        // point itself where they coincide).
        need(args, 4, "proj curve|surface x y z")?;
        let v = self.numbers(&args[1..4])?;
        let p = Point::new(v[0], v[1], v[2]);
        let found: Vec<(Point, String)> = match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(c)) => kentron::geometry::extrema::point_curve(p, &c)
                .into_iter()
                .enumerate()
                .map(|(i, q)| (q.point, format!(" parameter {} = {}", i + 1, format_real(q.parameter))))
                .collect(),
            Some(DrawObject::Surface(s)) => {
                let (u, w, _) = kentron::geometry::extrema::nearest_on_surface(p, &s, None);
                vec![(kentron::geometry::ParametricSurface::point(&s, u, w), format!(" Parameters: {} {}", format_real(u), format_real(w)))]
            }
            _ => return Err(DrawError::Failed(format!("{} is not a curve or a surface", args[0]))),
        };
        let mut out = String::new();
        for (i, (q, text)) in found.into_iter().enumerate() {
            let name = format!("ext_{}", i + 1);
            if p.distance(q) > tolerance::LINEAR {
                let line = Line::new(p, Direction::from_vector(q - p));
                self.set(&name, DrawObject::Curve(Curve::Trimmed(TrimmedCurve::new(Curve::Line(line), 0.0, p.distance(q)))));
            } else {
                self.set(&name, DrawObject::Point(q));
            }
            out += &format!("{name}{text}\n");
        }
        Ok(out)
    }
    /// DRAW command `project`: 47 uses in 46 converted scripts.
    ///
    /// Example: `project cx c s`
    pub fn project(&mut self, args: &[&str]) -> DrawResult {
        // `project result curve surface [-e tol] [-v ...] [-t tol] ...`: the curve projected
        // onto the surface (orthogonally, where it lies off it), in the surface's parameters
        // (the options tune OCCT's approximation).
        need(args, 3, "project result curve surface [options]")?;
        let curve = self.curve_var(args[1])?;
        let surface = self.surface_var(args[2])?;
        let projected = kentron::geometry::project::pcurve(&curve, &surface, 1e-7)
            .or_else(|| kentron::geometry::project::projected_pcurve(&curve, &surface, 1e-7))
            .ok_or_else(|| DrawError::Failed(format!("project: {} cannot be projected onto {}", args[1], args[2])))?;
        self.set(args[0], DrawObject::Curve2(Curve2::BSpline(projected)));
        Ok(String::new())
    }
    /// DRAW command `projonplane`: 1 uses in 1 converted scripts.
    ///
    /// Example: `projonplane result hyperbola plane`
    pub fn projonplane(&mut self, args: &[&str]) -> DrawResult {
        // `projonplane name curve plane [dx dy dz] [0|1]`: the curve projected onto the plane
        // along the direction (else its normal), with the curve's parametrization unless 0.
        need(args, 3, "projonplane name curve plane [dx dy dz] [0|1]")?;
        let curve = self.curve_var(args[1])?;
        let Surface::Plane(plane) = self.surface_var(args[2])? else {
            return Err(DrawError::Failed(" The surface must be a plane".into()));
        };
        let keep = match args.len() {
            4 => self.number(args[3])? != 0.0,
            7 => self.number(args[6])? != 0.0,
            _ => true,
        };
        let direction = if args.len() >= 6 {
            let v = self.numbers(&args[3..6])?;
            Direction::try_from_vector(Vector::new(v[0], v[1], v[2])).ok_or_else(|| DrawError::Failed("projonplane: null direction".into()))?
        } else {
            plane.normal()
        };
        let projected = kentron::geometry::project::curve_on_plane(&curve, &plane, direction, keep)
            .ok_or_else(|| DrawError::Failed("projonplane: the direction is parallel to the plane".into()))?;
        self.set(args[0], DrawObject::Curve(projected));
        Ok(String::new())
    }
    /// DRAW command `proximity`: 1 uses in 1 converted scripts.
    ///
    /// Example: `proximity s c -tol 0.01 -profile`
    pub fn proximity(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("proximity", args)
    }
    /// DRAW command `pscale`: 8 uses in 8 converted scripts.
    ///
    /// Example: `pscale result 0 0 0 1.2`
    pub fn pscale(&mut self, args: &[&str]) -> DrawResult {
        // `pscale name ... x y z factor`.
        need(args, 5, "pscale name ... x y z factor")?;
        let n = args.len();
        let v = self.numbers(&args[n - 4..])?;
        self.transform_vars(&args[..n - 4], &Transform::scaling(Point::new(v[0], v[1], v[2]), v[3]))
    }
    /// DRAW command `psphere`: 270 uses in 211 converted scripts.
    ///
    /// Example: `psphere s 15 -60 60 90`
    pub fn psphere(&mut self, args: &[&str]) -> DrawResult {
        // `psphere name [plane] radius [angle]` or `[latitude1 latitude2 [angle]]` (degrees).
        need(args, 2, "psphere name [plane] radius [angles]")?;
        let (frame, values) = self.primitive_frame(&args[1..])?;
        let sphere = match values[..] {
            [r] => Sphere::new(r),
            [r, angle] => Sphere::new(r).with_angle(angle.to_radians()),
            [r, a1, a2] => Sphere::new(r).with_latitudes(a1.to_radians(), a2.to_radians()),
            [r, a1, a2, angle] => {
                Sphere::new(r).with_latitudes(a1.to_radians(), a2.to_radians()).with_angle(angle.to_radians())
            }
            _ => return Err(DrawError::Failed("psphere: syntax error".into())),
        };
        let shape = sphere.with_frame(frame).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `ptorus`: 75 uses in 60 converted scripts.
    ///
    /// Example: `ptorus s 10 10 0 45 90`
    pub fn ptorus(&mut self, args: &[&str]) -> DrawResult {
        // `ptorus name [plane] R1 R2 [angle1 angle2] [angle]` (angles in degrees).
        need(args, 3, "ptorus name [plane] R1 R2 [angle1 angle2] [angle]")?;
        let (frame, values) = self.primitive_frame(&args[1..])?;
        let torus = match values[..] {
            [r1, r2] => Torus::new(r1, r2),
            [r1, r2, angle] => Torus::new(r1, r2).with_angle(angle.to_radians()),
            [r1, r2, a1, a2] => Torus::new(r1, r2).with_tube_angles(a1.to_radians(), a2.to_radians()),
            [r1, r2, a1, a2, angle] => {
                Torus::new(r1, r2).with_tube_angles(a1.to_radians(), a2.to_radians()).with_angle(angle.to_radians())
            }
            _ => return Err(DrawError::Failed("ptorus: syntax error".into())),
        };
        let shape = torus.with_frame(frame).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `pu`: 10 uses in 5 converted scripts.
    ///
    /// Example: `pu`
    pub fn pu(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("pu", args)
    }
    /// DRAW command `qcircle`: 5 uses in 3 converted scripts.
    ///
    /// Example: `qcircle c1 0 0 50`
    pub fn qcircle(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("qcircle", args)
    }
    /// DRAW command `qline`: 2 uses in 2 converted scripts.
    ///
    /// Example: `qline l1 -20 0 10 20`
    pub fn qline(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("qline", args)
    }
    /// DRAW command `reflectlines`: 5 uses in 5 converted scripts.
    ///
    /// Example: `reflectlines result a 0 1 0`
    pub fn reflectlines(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("reflectlines", args)
    }
    /// DRAW command `remknot`: 1 uses in 1 converted scripts.
    ///
    /// Example: `remknot result 2 0 0.01`
    pub fn remknot(&mut self, args: &[&str]) -> DrawResult {
        // `remknot curve index [mult [tol]]`: the B-spline curve's knot (from 1) lowered to the
        // multiplicity (0: removed) where the curve stays within the tolerance.
        need(args, 2, "remknot curve index [mult [tol]]")?;
        let index = (self.number(args[1])? as usize).saturating_sub(1);
        let mult = match args.get(2) { Some(m) => self.number(m)? as usize, None => 0 };
        let tol = match args.get(3) { Some(t) => self.number(t)?, None => 1e-7 };
        let done = match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(Curve::BSpline(mut b))) => { let ok = b.remove_knot(index, mult, tol); self.set(args[0], DrawObject::Curve(Curve::BSpline(b))); ok }
            Some(DrawObject::Curve2(Curve2::BSpline(mut b))) => { let ok = b.remove_knot(index, mult, tol); self.set(args[0], DrawObject::Curve2(Curve2::BSpline(b))); ok }
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline curve", args[0]))),
        };
        Ok(if done { String::new() } else { "Remove knots failed".into() })
    }
    /// DRAW command `removeinternals`: 6 uses in 4 converted scripts.
    ///
    /// Example: `removeinternals r`
    pub fn removeinternals(&mut self, args: &[&str]) -> DrawResult {
        // `removeinternals shape [force]`: the shape without its sub-shapes of internal
        // orientation.
        need(args, 1, "removeinternals shape [force]")?;
        let shape = self.shape_var(args[0])?;
        let internal: Vec<Shape> = self.model.explore_all(shape).filter(|s| s.orientation == Orientation::Internal).collect();
        if internal.is_empty() {
            return Ok(String::new());
        }
        let mut sub = kentron::topology::Substitution::new();
        for s in internal {
            sub.remove(s);
        }
        let result = sub.apply(&mut self.model, shape).ok_or_else(|| DrawError::Failed("removeinternals: nothing left".into()))?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `removeloc`: 1 uses in 1 converted scripts.
    ///
    /// Example: `removeloc r a …`
    pub fn removeloc(&mut self, args: &[&str]) -> DrawResult {
        // `removeloc result shape`: a copy with its placement put into its geometry.
        need(args, 2, "removeloc result shape")?;
        let shape = self.shape_var(args[1])?;
        let copy = modify::copy(&mut self.model, shape, false)?;
        self.set_shape_var(args[0], copy)
    }
    /// DRAW command `renamevar`: 122 uses in 113 converted scripts.
    ///
    /// Example: `renamevar result_1 result`
    pub fn renamevar(&mut self, args: &[&str]) -> DrawResult {
        // `renamevar from to [from to ...]`: variables renamed.
        if args.is_empty() || args.len() % 2 != 0 {
            return Err(DrawError::Failed("renamevar from to [from to ...]".into()));
        }
        for pair in args.chunks(2) {
            let value = self.take_var(pair[0])?;
            self.set(pair[1], value);
        }
        Ok(String::new())
    }
    /// DRAW command `repeat`: 54 uses in 49 converted scripts.
    ///
    /// Example: `repeat 199 plane p_$i 0 $i*100 0 0 1 0 mkface f_$i p_$i set bsres [bsection s_$i a f_$i] i`
    pub fn repeat(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("repeat", args)
    }
    /// DRAW command `restore`: 8680 uses in 5956 converted scripts.
    ///
    /// Example: `restore … s`
    pub fn restore(&mut self, args: &[&str]) -> DrawResult {
        // `restore file name`: a shape, or a point, curve or surface DRAW saved.
        let [file, name] = args else { return self.not_implemented("restore", args) };
        let bytes = std::fs::read(file).map_err(|e| DrawError::Failed(format!("restore: {file}: {e}")))?;
        use kentron::io::brep::SavedGeometry;
        let object = match kentron::io::brep::read_saved_geometry(&mut bytes.as_slice())? {
            Some(SavedGeometry::Point(p)) => DrawObject::Point(p),
            Some(SavedGeometry::Point2(p)) => DrawObject::Point2(p),
            Some(SavedGeometry::Curve(c)) => DrawObject::Curve(c),
            Some(SavedGeometry::Curve2(c)) => DrawObject::Curve2(c),
            Some(SavedGeometry::Surface(s)) => DrawObject::Surface(s),
            None => {
                let shape = kentron::io::brep::read(&mut self.model, &mut bytes.as_slice())?;
                return self.set_shape_var(name, shape);
            }
        };
        self.set(name, object);
        Ok(String::new())
    }
    /// DRAW command `save`.
    ///
    /// Example: `save result result.brep`
    pub fn save(&mut self, args: &[&str]) -> DrawResult {
        // `save name file`: the shape written as a BRep file (as `restore` reads it).
        let [name, file] = args else { return self.not_implemented("save", args) };
        let shape = self.shape_var(name)?;
        let mut out = std::fs::File::create(file).map_err(|e| DrawError::Failed(format!("save: {file}: {e}")))?;
        kentron::io::brep::write(&self.model, shape, &mut out, &Default::default())?;
        Ok(format!("{name} saved\n"))
    }
    /// DRAW command `result`: 8 uses in 5 converted scripts.
    ///
    /// Example: `result`
    pub fn result(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("result", args)
    }
    /// DRAW command `reverse`: 16 uses in 16 converted scripts.
    ///
    /// Example: `reverse result`
    pub fn reverse(&mut self, args: &[&str]) -> DrawResult {
        // `reverse name ...`: reverses the curves' parametrization (a shape's orientation).
        for name in args {
            let object = match self.get(name).cloned() {
                Some(DrawObject::Curve(mut c)) => { c.reverse(); DrawObject::Curve(c) }
                Some(DrawObject::Curve2(mut c)) => { c.reverse(); DrawObject::Curve2(c) }
                Some(DrawObject::Shape(s)) => DrawObject::Shape(s.reversed()),
                _ => return self.not_implemented("reverse", args),
            };
            self.set(name, object);
        }
        Ok(String::new())
    }
    /// DRAW command `revol`: 85 uses in 79 converted scripts.
    ///
    /// Example: `revol s w 0 0 0 0 0 1 270`
    pub fn revol(&mut self, args: &[&str]) -> DrawResult {
        // `revol result shape px py pz dx dy dz angle [copy]`: the angle is in degrees. The copy
        // flag changes nothing in kentron (the swept geometry is always new).
        need(args, 9, "revol result shape px py pz dx dy dz angle [copy]")?;
        let profile = self.shape_var(args[1])?;
        let values = args[2..9].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let direction = Direction::try_from_vector(Vector::new(values[3], values[4], values[5]))
            .ok_or_else(|| DrawError::Failed("revol: null axis direction".into()))?;
        let axis = Axis::new(Point::new(values[0], values[1], values[2]), direction);
        let result = sweep::revolve(&mut self.model, profile, axis, values[6].to_radians())?;
        self.set_shape_var(args[0], result.shape)
    }
    /// DRAW command `revsurf`: 30 uses in 30 converted scripts.
    ///
    /// Example: `revsurf ss cc 0 0 0 0 0 1`
    pub fn revsurf(&mut self, args: &[&str]) -> DrawResult {
        // `revsurf name curve x y z dx dy dz`: the surface swept by the curve (or an edge's)
        // around an axis.
        need(args, 8, "revsurf name curve x y z dx dy dz")?;
        let curve = self.curve_or_edge(args[1])?;
        let v = self.numbers(&args[2..8])?;
        let d = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("revsurf: null direction".into()))?;
        let axis = Axis::new(Point::new(v[0], v[1], v[2]), d);
        self.set(args[0], DrawObject::Surface(Surface::Revolution(RevolutionSurface::new(curve, axis))));
        Ok(String::new())
    }
    /// DRAW command `rotate`: 2 uses in 1 converted scripts.
    ///
    /// Example: `rotate ee2 0 0 100 0 1 0 10`
    pub fn rotate(&mut self, args: &[&str]) -> DrawResult {
        // `rotate name ... x y z dx dy dz angle` (the angle in degrees).
        need(args, 8, "rotate name ... x y z dx dy dz angle")?;
        let n = args.len();
        let v = self.numbers(&args[n - 7..])?;
        let d = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("rotate: null direction".into()))?;
        let t = Transform::rotation(Axis::new(Point::new(v[0], v[1], v[2]), d), v[6].to_radians());
        self.transform_vars(&args[..n - 7], &t)
    }
    /// DRAW command `sameparameter`: 7 uses in 5 converted scripts.
    ///
    /// Example: `sameparameter result`
    pub fn sameparameter(&mut self, args: &[&str]) -> DrawResult {
        // `sameparameter [result] shape [tolerance]`: kentron keeps edges same-parameter as it
        // builds them; this checks the shape is (fails if not) and names it `result`. (A
        // positive last argument is the tolerance, as OCCT reads it.)
        need(args, 1, "sameparameter [result] shape [tolerance]")?;
        let tolerance = args.len() >= 2 && args[args.len() - 1].parse::<f64>().is_ok_and(|t| t > 0.0);
        let input = if tolerance { args[args.len() - 2] } else { args[args.len() - 1] };
        let shape = self.shape_var(input)?;
        let report = validity::check(&self.model, shape);
        let problems: Vec<String> = report.errors().map(|p| p.description.clone()).filter(|d| d.contains("pcurve")).collect();
        if !problems.is_empty() {
            return self.not_implemented("sameparameter", args);
        }
        if args[0] != input {
            self.set_shape_var(args[0], shape)?;
        }
        Ok(String::new())
    }
    /// DRAW command `samerange`: 1 uses in 1 converted scripts.
    ///
    /// Example: `samerange comp_2`
    pub fn samerange(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("samerange", args)
    }
    /// DRAW command `savehistory`: 47 uses in 43 converted scripts.
    ///
    /// Example: `savehistory hist`
    /// The history saved as `name`.
    fn draw_history(&self, name: &str) -> Result<DrawHistory, DrawError> {
        match self.get(name) {
            Some(DrawObject::History(h)) => Ok(h.clone()),
            _ => Err(DrawError::Failed(format!("{name} is not a history"))),
        }
    }
    pub fn savehistory(&mut self, args: &[&str]) -> DrawResult {
        // `savehistory name`: the last operation's history (of those keeping one).
        need(args, 1, "savehistory name")?;
        self.set(args[0], DrawObject::History(self.last_history.clone()));
        Ok(String::new())
    }
    /// DRAW command `scaleshape`: 2 uses in 2 converted scripts.
    ///
    /// Example: `scaleshape result a 10`
    pub fn scaleshape(&mut self, args: &[&str]) -> DrawResult {
        // `scaleshape result shape factor`: scaled about the origin.
        need(args, 3, "scaleshape result shape factor")?;
        let shape = self.shape_var(args[1])?;
        let t = Transform::scaling(Point::ORIGIN, self.number(args[2])?);
        let scaled = modify::transform(&mut self.model, shape, &t, Default::default())?;
        self.set_shape_var(args[0], scaled)
    }
    /// DRAW command `scalexyz`: 5 uses in 4 converted scripts.
    ///
    /// Example: `scalexyz result a 2 3 4`
    pub fn scalexyz(&mut self, args: &[&str]) -> DrawResult {
        // `scalexyz result shape kx ky kz`: scaled along the axes by each factor.
        need(args, 5, "scalexyz result shape kx ky kz")?;
        let shape = self.shape_var(args[1])?;
        let k = self.numbers(&args[2..5])?;
        let m = kentron::math::Matrix3::from_rows([[k[0], 0.0, 0.0], [0.0, k[1], 0.0], [0.0, 0.0, k[2]]]);
        let t = kentron::math::AffineTransform::new(m, Vector::ZERO);
        let scaled = modify::transform_affine(&mut self.model, shape, &t)?;
        self.set_shape_var(args[0], scaled)
    }
    /// DRAW command `segment`: 2 uses in 2 converted scripts.
    ///
    /// Example: `segment result -4.76283758432072e-05 1.06458532027221`
    pub fn segment(&mut self, args: &[&str]) -> DrawResult {
        // `segment curve u1 u2`: the (B-spline or Bézier) curve restricted to [u1, u2]; other
        // curves trimmed.
        need(args, 3, "segment curve u1 u2")?;
        let (u1, u2) = (self.number(args[1])?, self.number(args[2])?);
        match self.get(args[0]).cloned() {
            Some(DrawObject::Curve(Curve::BSpline(mut b))) => {
                b.segment(u1, u2);
                self.set(args[0], DrawObject::Curve(Curve::BSpline(b)));
            }
            Some(DrawObject::Curve(c)) => self.set(args[0], DrawObject::Curve(Curve::Trimmed(TrimmedCurve::new(c, u1, u2)))),
            _ => return self.not_implemented("segment", args),
        }
        Ok(String::new())
    }
    /// DRAW command `segsur`: 1 uses in 1 converted scripts.
    ///
    /// Example: `segsur s 0.0, 1.0, 0.49999999999999988898, 0.66666666666666662966`
    pub fn segsur(&mut self, args: &[&str]) -> DrawResult {
        // `segsur surface u1 u2 v1 v2 [utol [vtol]]`: the B-spline or Bézier surface restricted to
        // [u1, u2] x [v1, v2].
        need(args, 5, "segsur surface u1 u2 v1 v2")?;
        let v = self.numbers(&args[1..5])?;
        let s = match self.surface_var(args[0])? {
            Surface::BSpline(mut b) => { b.segment(v[0], v[1], v[2], v[3]); Surface::BSpline(b) }
            Surface::Bezier(mut b) => { b.segment(v[0], v[1], v[2], v[3]); Surface::Bezier(b) }
            _ => return Err(DrawError::Failed(format!("{} is not a B-spline or Bézier surface", args[0]))),
        };
        self.set(args[0], DrawObject::Surface(s));
        Ok(String::new())
    }
    /// DRAW command `setflags`: 79 uses in 41 converted scripts.
    ///
    /// Example: `setflags b locked`
    pub fn setflags(&mut self, args: &[&str]) -> DrawResult {
        // `setflags shape [-]flag...`: OCCT's locked/modified/checked/... marks, which guard
        // against or record in-place changes. kentron's operations do not change their
        // arguments, so the flags have nothing to do.
        let [shape, ..] = args else { return self.not_implemented("setflags", args) };
        self.shape_var(shape)?;
        Ok(String::new())
    }
    /// DRAW command `setsweep`: 31 uses in 31 converted scripts.
    ///
    /// Example: `setsweep -FR`
    pub fn setsweep(&mut self, args: &[&str]) -> DrawResult {
        // `setsweep -CF|-DT|-FR|...`: the trihedron's mode. The rotation-minimizing frame
        // kentron sweeps with is OCCT's corrected Frenet (`-CF`) and what its discrete
        // trihedron (`-DT`) approximates; Frenet's (`-FR`) is the same along a plane curve.
        // Fixed, constant-binormal, guided and surface modes are not supported.
        match args {
            ["-CF"] | ["-DT"] | ["-FR"] => {
                self.sweep_mode = args.iter().map(|a| a.to_string()).collect();
                Ok(String::new())
            }
            _ => self.not_implemented("setsweep", args),
        }
    }
    /// DRAW command `settolerance`: 145 uses in 45 converted scripts.
    ///
    /// Example: `settolerance a 0.1`
    pub fn settolerance(&mut self, args: &[&str]) -> DrawResult {
        // `settolerance shape [v|e|f|a] value` or `settolerance shape [v|e|f|a] min max`: sets
        // the tolerances of the sub-shapes, or brings them into [min, max].
        need(args, 2, "settolerance shape [v|e|f|a] value | min max")?;
        let shape = self.shape_var(args[0])?;
        let (target, rest) = match args[1].to_ascii_lowercase().as_str() {
            "v" => (modify::ToleranceTarget::Vertices, &args[2..]),
            "e" => (modify::ToleranceTarget::Edges, &args[2..]),
            "f" => (modify::ToleranceTarget::Faces, &args[2..]),
            "a" => (modify::ToleranceTarget::All, &args[2..]),
            _ => (modify::ToleranceTarget::All, &args[1..]),
        };
        match self.numbers(rest)?[..] {
            [value] => modify::set_tolerances(&mut self.model, shape, value, target)?,
            [min, max] => {
                modify::limit_tolerances(&mut self.model, shape, min, max, target)?;
            }
            _ => return self.not_implemented("settolerance", args),
        }
        Ok(String::new())
    }
    /// DRAW command `setunotperiodic`: 1 uses in 1 converted scripts.
    ///
    /// Example: `setunotperiodic bs`
    pub fn setunotperiodic(&mut self, args: &[&str]) -> DrawResult {
        // `setunotperiodic surface ...`: the B-spline surfaces made non-periodic in u.
        need(args, 1, "setunotperiodic surface ...")?;
        for name in args {
            let s = match self.surface_var(name)? {
                Surface::BSpline(mut b) => { b.set_u_not_periodic(); Surface::BSpline(b) }
                other => other,
            };
            self.set(name, DrawObject::Surface(s));
        }
        Ok(String::new())
    }
    /// DRAW command `sewing`: 752 uses in 747 converted scripts.
    ///
    /// Example: `sewing se b1_1 b1_2 b1_3 b1_4 b1_5 res_1 b2_5`
    pub fn sewing(&mut self, args: &[&str]) -> DrawResult {
        // `sewing result [tolerance | +t tolerance] shapes... [+mint t] [+maxt t] [switches]`: of the
        // switches, `-c` turns cutting off and `-n` asks for non-manifold sewing; the others
        // (analysis, face mode, floating edges, same-parameter processing) have no counterpart.
        let Some((&result, rest)) = args.split_first() else { return self.not_implemented("sewing", args) };
        let mut options = kentron::sew::SewOptions::default();
        let mut shapes = Vec::new();
        let mut i = 0;
        while i < rest.len() {
            match rest[i] {
                "+t" => {
                    i += 1;
                    options.tolerance = rest.get(i).and_then(|t| t.parse().ok()).ok_or_else(|| DrawError::Failed("sewing: +t needs a tolerance".into()))?;
                }
                "+mint" | "+maxt" => i += 1,
                a if a.starts_with('+') => {}
                "-c" => options.cutting = false,
                "-n" => options.non_manifold = true,
                a if a.starts_with('-') && a.parse::<f64>().is_err() => {}
                a => match a.parse::<f64>() {
                    Ok(t) if i == 0 => options.tolerance = t,
                    Ok(_) => {} // the minimal and maximal tolerances
                    Err(_) => shapes.push(self.shape_var(a)?),
                },
            }
            i += 1;
        }
        let sewn = kentron::sew::sew(&mut self.model, &shapes, &options)?;
        self.set_shape_var(result, sewn)
    }
    /// DRAW command `shape`: 40 uses in 32 converted scripts.
    ///
    /// Example: `shape s So`
    pub fn shape(&mut self, args: &[&str]) -> DrawResult {
        // `shape name Sh|So|CS|C|W|F`: an empty shape of the type.
        need(args, 2, "shape name type")?;
        let kind = match args[1].to_ascii_uppercase().as_str() {
            "C" => NodeKind::Compound,
            "CS" => NodeKind::CompSolid,
            "SO" => NodeKind::Solid,
            "SH" => NodeKind::Shell,
            "F" => NodeKind::Face,
            "W" => NodeKind::Wire,
            _ => return self.not_implemented("shape", args),
        };
        let shape = self.model.new_empty(kind);
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `side`: 6 uses in 3 converted scripts.
    ///
    /// Example: `side left`
    pub fn side(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("side", args)
    }
    /// DRAW command `smirror`: 12 uses in 12 converted scripts.
    ///
    /// Example: `smirror result 0 0 0 1 0 0`
    pub fn smirror(&mut self, args: &[&str]) -> DrawResult {
        // `smirror name ... x y z nx ny nz`: the symmetry about a plane.
        need(args, 7, "smirror name ... x y z nx ny nz")?;
        let n = args.len();
        let v = self.numbers(&args[n - 6..])?;
        let d = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("smirror: null direction".into()))?;
        self.transform_vars(&args[..n - 6], &Transform::mirror_plane(Point::new(v[0], v[1], v[2]), d))
    }
    /// DRAW command `solidorientation`: 1 uses in 1 converted scripts.
    ///
    /// Example: `solidorientation result`
    pub fn solidorientation(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("solidorientation", args)
    }
    /// DRAW command `sphere`: 35 uses in 27 converted scripts.
    ///
    /// Example: `sphere sph_f1 0 0 0 0 0 1 550`
    pub fn sphere(&mut self, args: &[&str]) -> DrawResult {
        // `sphere name [x y z [nx ny nz [ux uy uz]]] radius`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(1));
        let frame = if frame_args.is_empty() { Frame::WORLD } else { Draw::curve_frame(frame_args)? };
        self.set(args[0], DrawObject::Surface(Surface::Sphere(kentron::geometry::surface::Sphere::new(frame, r[0]))));
        Ok(String::new())
    }
    /// DRAW command `splitc12d`: 1 uses in 1 converted scripts.
    ///
    /// Example: `splitc12d pc 0`
    pub fn splitc12d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("splitc12d", args)
    }
    /// DRAW command `splitshape`: 13 uses in 13 converted scripts.
    ///
    /// Example: `splitshape result f f e`
    pub fn splitshape(&mut self, args: &[&str]) -> DrawResult {
        // `splitshape result shape [face edge|wire|compound ...]... | [edges ...]`: the faces
        // divided along edges lying on them (on the face named before them; without one, on
        // whichever faces they lie on). `@` pairs (edges on edges) are not supported. The
        // result is made from a copy.
        need(args, 3, "splitshape result shape [face edges ...]... | edges ...")?;
        if args.iter().any(|a| a.starts_with('@') || a.starts_with('-') || a.starts_with('.')) {
            return self.not_implemented("splitshape", args);
        }
        let shape = self.shape_var(args[1])?;
        let copy = kentron::modify::copy(&mut self.model, shape, false)?;
        // (A face named: its copy, the face of the copy with its area and centre.)
        let copied_faces = self.model.unique_subshapes(copy, NodeKind::Face);
        let props = |m: &kentron::topology::Model, f: Shape| {
            let p = kentron::properties::area(m, f);
            (p.mass, p.center_of_mass)
        };
        // Groups: a face and the edges after it, or edges with no face.
        let mut groups: Vec<(Option<Shape>, Vec<Shape>)> = Vec::new();
        for name in &args[2..] {
            let s = self.shape_var(name)?;
            if s.kind() == NodeKind::Face {
                let (area, center) = props(&self.model, s);
                let face = copied_faces
                    .iter()
                    .copied()
                    .find(|&f| {
                        let (a, c) = props(&self.model, f);
                        (a - area).abs() <= 1e-9 * (1.0 + area) && c.distance(center) <= 1e-7 * (1.0 + area.sqrt())
                    })
                    .ok_or_else(|| DrawError::Failed(format!("{name} is not a face of {}", args[1])))?;
                groups.push((Some(face), Vec::new()));
            } else {
                let edges = self.model.unique_subshapes(s, NodeKind::Edge);
                match groups.last_mut() {
                    Some((_, list)) => list.extend(edges),
                    None => groups.push((None, edges)),
                }
            }
        }
        let mut result = copy;
        for (face, edges) in groups {
            result = kentron::divide::divide_by_edges(&mut self.model, result, face, &edges)?;
        }
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `ssolid`: 4 uses in 1 converted scripts.
    ///
    /// Example: `ssolid s1 s1`
    pub fn ssolid(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("ssolid", args)
    }
    /// DRAW command `subshape`: 2 uses in 2 converted scripts.
    ///
    /// Example: `subshape sh F 379`
    pub fn subshape(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("subshape", args)
    }
    /// DRAW command `svalue`: 7 uses in 3 converted scripts.
    ///
    /// Example: `svalue otobz 0.2 0.2 x y z dux duy duz dvx dvy dvz`
    pub fn svalue(&mut self, args: &[&str]) -> DrawResult {
        // `svalue surface u v x y z [dux duy duz dvx dvy dvz [d2u.. d2v.. duv..]]`: the point and
        // derivatives at (u, v), into numeric variables.
        need(args, 6, "svalue surface u v x y z [...]")?;
        use kentron::geometry::ParametricSurface;
        let s = self.surface_var(args[0])?;
        let (u, v) = (self.number(args[1])?, self.number(args[2])?);
        let d = s.d2(u, v);
        let values = [d.point.x, d.point.y, d.point.z, d.du.x, d.du.y, d.du.z, d.dv.x, d.dv.y, d.dv.z, d.duu.x, d.duu.y, d.duu.z, d.dvv.x, d.dvv.y, d.dvv.z, d.duv.x, d.duv.y, d.duv.z];
        for (name, value) in args[3..].iter().zip(values) {
            self.set(name, DrawObject::Number(value));
        }
        Ok(String::new())
    }
    /// DRAW command `tclean`: 181 uses in 167 converted scripts.
    ///
    /// Example: `tclean a`
    pub fn tclean(&mut self, args: &[&str]) -> DrawResult {
        // `tclean [-geom] [-force] name ...`: removes the meshes stored on the shapes.
        for name in args.iter().filter(|a| !a.starts_with('-')) {
            let shape = self.shape_var(name)?;
            self.model.clear_meshes(shape);
        }
        Ok(String::new())
    }
    /// DRAW command `tcopy`: 80 uses in 42 converted scripts.
    ///
    /// Example: `tcopy result_1 result`
    pub fn tcopy(&mut self, args: &[&str]) -> DrawResult {
        // `tcopy name copy`: a deep copy.
        need(args, 2, "tcopy name copy")?;
        let shape = self.shape_var(args[0])?;
        let copy = modify::copy(&mut self.model, shape, false)?;
        self.set_shape_var(args[1], copy)
    }
    /// DRAW command `tessellate`: 1 uses in 1 converted scripts.
    ///
    /// Example: `tessellate r result_1 50 50`
    pub fn tessellate(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("tessellate", args)
    }
    /// DRAW command `thickshell`: 53 uses in 53 converted scripts.
    ///
    /// Example: `thickshell result s 5`
    pub fn thickshell(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("thickshell", args)
    }
    /// DRAW command `thrusections`: 293 uses in 291 converted scripts.
    ///
    /// Example: `thrusections r 0 0 w_1 w_2`
    pub fn thrusections(&mut self, args: &[&str]) -> DrawResult {
        // `thrusections [-N] result issolid isruled sections... [-safe]` (`-N` takes the
        // sections as they are, without OCCT's compatibility checks; `-safe` skips its
        // in-place modifications, which changes nothing here). The history is kept.
        let as_given = args.contains(&"-N");
        let args: Vec<&str> = args.iter().copied().filter(|a| !a.starts_with('-')).collect();
        if args.len() < 5 {
            return self.not_implemented("thrusections", &args);
        }
        let flag = |a: &str| a.parse::<f64>().map(|x| x != 0.0).map_err(|_| DrawError::Failed(format!("thrusections: {a} is not 0 or 1")));
        let options = kentron::loft::LoftOptions { solid: flag(args[1])?, ruled: flag(args[2])?, matching: if as_given { kentron::loft::Matching::AsGiven } else { kentron::loft::Matching::ByRays }, ..Default::default() };
        let sections = args[3..].iter().map(|a| self.shape_var(a)).collect::<Result<Vec<_>, _>>()?;
        let result = kentron::loft::loft(&mut self.model, &sections, &options)?;
        let generated = sections.iter().flat_map(|&s| {
            let parts = self.model.unique_subshapes(s, NodeKind::Edge).into_iter().chain(self.model.unique_subshapes(s, NodeKind::Vertex));
            parts.map(|p| (p.node, result.generated(p))).filter(|(_, g)| !g.is_empty()).collect::<Vec<_>>()
        }).collect();
        self.last_history = DrawHistory { generated, ..Default::default() };
        // (OCCT flags its solids closed.)
        if result.shape.kind() == NodeKind::Solid {
            self.closed_solids.insert(result.shape.node);
        }
        self.set_shape_var(args[0], result.shape)
    }
    /// DRAW command `tmirror`: 5 uses in 5 converted scripts.
    ///
    /// Example: `tmirror b 0 0 0 1 0 0`
    pub fn tmirror(&mut self, args: &[&str]) -> DrawResult {
        // `tmirror name ... x y z nx ny nz [-copy] [-copymesh]`: mirrors the shapes in the
        // plane through the point with the normal (rewriting the geometry, as a mirror is no
        // placement).
        let options = transform_options(args);
        let args = &without_copy_flags(args);
        need(args, 7, "tmirror name ... x y z nx ny nz")?;
        let n = args.len();
        let v = self.numbers(&args[n - 6..])?;
        let normal = Direction::try_from_vector(Vector::new(v[3], v[4], v[5])).ok_or_else(|| DrawError::Failed("tmirror: null normal".into()))?;
        let t = Transform::mirror_plane(Point::new(v[0], v[1], v[2]), normal);
        for name in &args[..n - 6] {
            let shape = self.shape_var(name)?;
            let mirrored = modify::transform(&mut self.model, shape, &t, options)?;
            self.set_shape_var(name, mirrored)?;
        }
        Ok(String::new())
    }
    /// DRAW command `to3d`: 24 uses in 4 converted scripts.
    ///
    /// Example: `to3d ee1 e1`
    pub fn to3d(&mut self, args: &[&str]) -> DrawResult {
        // `to3d name curve2d [plane]`: the 2D curve placed on the plane (XY if none).
        need(args, 2, "to3d name curve2d [plane]")?;
        let c = match self.get(args[1]) {
            Some(DrawObject::Curve2(c)) => c.clone(),
            _ => return self.not_implemented("to3d", args),
        };
        let frame = match args.get(2) {
            Some(p) => match self.surface_var(p)? {
                Surface::Plane(pl) => *pl.frame(),
                _ => return Err(DrawError::Failed(format!("{p} is not a plane"))),
            },
            None => Frame::WORLD,
        };
        self.set(args[0], DrawObject::Curve(c.to_3d(&frame)));
        Ok(String::new())
    }
    /// DRAW command `tocanon`: 2 uses in 2 converted scripts.
    ///
    /// Example: `tocanon r c`
    pub fn tocanon(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("tocanon", args)
    }
    /// DRAW command `tolerance`: 23 uses in 16 converted scripts.
    ///
    /// Example: `tolerance a`
    pub fn tolerance(&mut self, args: &[&str]) -> DrawResult {
        // `tolerance shape`: the largest, mean and smallest tolerances of its sub-shapes, all
        // together, then its faces', edges' and vertices' (each counted per use, as OCCT's).
        need(args, 1, "tolerance shape")?;
        let shape = self.shape_var(args[0])?;
        let tolerances = |kinds: &[NodeKind]| -> Vec<f64> {
            kinds.iter().flat_map(|&k| self.model.explore(shape, k).map(|s| self.model.tolerance(s)).collect::<Vec<_>>()).collect()
        };
        let line = |values: &[f64]| {
            if values.is_empty() {
                return "MAX=0 AVG=0 MIN=0".to_string();
            }
            let max = values.iter().copied().fold(f64::MIN, f64::max);
            let min = values.iter().copied().fold(f64::MAX, f64::min);
            format!("MAX={max} AVG={} MIN={min}", values.iter().sum::<f64>() / values.len() as f64)
        };
        let mut out = format!("Tolerance {}\n", line(&tolerances(&[NodeKind::Face, NodeKind::Edge, NodeKind::Vertex])));
        for (name, kind) in [("FACE    ", NodeKind::Face), ("EDGE    ", NodeKind::Edge), ("VERTEX  ", NodeKind::Vertex)] {
            out += &format!("{name}: {}\n", line(&tolerances(&[kind])));
        }
        Ok(out)
    }
    /// DRAW command `topoload`: 7 uses in 7 converted scripts.
    ///
    /// Example: `topoload res`
    pub fn topoload(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("topoload", args)
    }
    /// DRAW command `torus`: 5 uses in 5 converted scripts.
    ///
    /// Example: `torus a 10 10`
    pub fn torus(&mut self, args: &[&str]) -> DrawResult {
        // `torus name [x y z [nx ny nz [ux uy uz]]] major minor`.
        let v = self.numbers(&args[1..])?;
        let (frame_args, r) = v.split_at(v.len().saturating_sub(2));
        let frame = if frame_args.is_empty() { Frame::WORLD } else { Draw::curve_frame(frame_args)? };
        self.set(args[0], DrawObject::Surface(Surface::Torus(kentron::geometry::surface::Torus::new(frame, r[0], r[1]))));
        Ok(String::new())
    }
    /// DRAW command `translate`: 22 uses in 21 converted scripts.
    ///
    /// Example: `translate ee2 0 0 100`
    pub fn translate(&mut self, args: &[&str]) -> DrawResult {
        // `translate name ... dx dy dz`.
        need(args, 4, "translate name ... dx dy dz")?;
        let n = args.len();
        let v = self.numbers(&args[n - 3..])?;
        self.transform_vars(&args[..n - 3], &Transform::translation(Vector::new(v[0], v[1], v[2])))
    }
    /// DRAW command `treverse`: 2 uses in 2 converted scripts.
    ///
    /// Example: `treverse s`
    pub fn treverse(&mut self, args: &[&str]) -> DrawResult {
        // `treverse name ...`: reverses the shapes' orientations.
        for name in args {
            let shape = self.shape_var(name)?;
            self.set_shape_var(name, shape.reversed())?;
        }
        Ok(String::new())
    }
    /// DRAW command `triangles`: 50 uses in 50 converted scripts.
    ///
    /// Example: `triangles result`
    pub fn triangles(&mut self, args: &[&str]) -> DrawResult {
        // `triangles shape`: shows the shape's triangulation in DRAW's viewer; nothing to
        // compute.
        need(args, 1, "triangles shape")?;
        Ok(String::new())
    }
    /// DRAW command `tricheck`: 8 uses in 5 converted scripts.
    ///
    /// Example: `tricheck r_1`
    pub fn tricheck(&mut self, args: &[&str]) -> DrawResult {
        // `tricheck shape`: problems of the faces' meshes (nothing if none): triangles with
        // nodes out of range or repeated, mesh edges shared by more than two triangles,
        // nodes in no triangle.
        need(args, 1, "tricheck shape")?;
        let shape = self.shape_var(args[0])?;
        let mut out = String::new();
        for (index, face) in self.model.unique_subshapes(shape, NodeKind::Face).into_iter().enumerate() {
            let Some(mesh) = self.model.face_mesh(face) else { continue };
            let n = mesh.nodes().len();
            let mut used = vec![false; n];
            let mut links: std::collections::HashMap<(usize, usize), usize> = std::collections::HashMap::new();
            for (t, tri) in mesh.triangles().iter().enumerate() {
                let tri = tri.map(|i| i as usize);
                if tri.iter().any(|&i| i >= n) || tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                    out += &format!("face {} triangle {}: bad nodes {:?}\n", index + 1, t + 1, tri);
                    continue;
                }
                for k in 0..3 {
                    used[tri[k]] = true;
                    let (a, b) = (tri[k], tri[(k + 1) % 3]);
                    *links.entry((a.min(b), a.max(b))).or_default() += 1;
                }
            }
            let shared = links.values().filter(|&&c| c > 2).count();
            if shared > 0 {
                out += &format!("face {}: {shared} links shared by more than two triangles\n", index + 1);
            }
            let free = used.iter().filter(|u| !**u).count();
            if free > 0 {
                out += &format!("face {}: Free_nodes {free}\n", index + 1);
            }
        }
        Ok(out)
    }
    /// DRAW command `trim`: 179 uses in 105 converted scripts.
    ///
    /// Example: `trim cc1 c1 0 3.14`
    pub fn trim(&mut self, args: &[&str]) -> DrawResult {
        // `trim name basis [u1 u2 [v1 v2]]`: a trimmed curve or surface (of the basis of an
        // already trimmed one); without a range, the untrimmed basis.
        need(args, 2, "trim name basis [u1 u2 [v1 v2]]")?;
        let v = self.numbers(&args[2..])?;
        let object = match (self.get(args[1]).cloned(), v.len()) {
            (Some(DrawObject::Curve(c)), n @ (0 | 2)) => {
                let basis = match c { Curve::Trimmed(t) => t.basis().clone(), c => c };
                DrawObject::Curve(if n == 0 { basis } else { Curve::Trimmed(TrimmedCurve::new(basis, v[0], v[1])) })
            }
            (Some(DrawObject::Curve2(c)), n @ (0 | 2)) => {
                let basis = match c { Curve2::Trimmed(t) => t.basis().clone(), c => c };
                DrawObject::Curve2(if n == 0 { basis } else { Curve2::Trimmed(TrimmedCurve2::new(basis, v[0], v[1])) })
            }
            (Some(DrawObject::Surface(s)), n @ (0 | 4)) => {
                let basis = match s { Surface::Trimmed(t) => t.basis().clone(), s => s };
                DrawObject::Surface(if n == 0 { basis } else { Surface::Trimmed(TrimmedSurface::new(basis, v[0], v[1], v[2], v[3])) })
            }
            _ => return self.not_implemented("trim", args),
        };
        self.set(args[0], object);
        Ok(String::new())
    }
    /// DRAW command `trimu`: 1 uses in 1 converted scripts.
    ///
    /// Example: `trimu ts ss 0 pi`
    pub fn trimu(&mut self, args: &[&str]) -> DrawResult {
        // `trimu name surface u1 u2`: the surface trimmed in u only.
        need(args, 4, "trimu name surface u1 u2")?;
        // (Of an already trimmed surface, its basis keeping the other direction's trim, as
        // OCCT's.)
        let surface = self.surface_var(args[1])?;
        let (u1, u2) = (self.number(args[2])?, self.number(args[3])?);
        let (_, _, v1, v2) = kentron::geometry::ParametricSurface::bounds(&surface);
        let surface = match surface {
            Surface::Trimmed(t) => t.basis().clone(),
            s => s,
        };
        self.set(args[0], DrawObject::Surface(Surface::Trimmed(TrimmedSurface::new(surface, u1, u2, v1, v2))));
        Ok(String::new())
    }
    /// DRAW command `trimv`: 25 uses in 12 converted scripts.
    ///
    /// Example: `trimv c1 c1 0 50`
    pub fn trimv(&mut self, args: &[&str]) -> DrawResult {
        // `trimv name surface v1 v2`: the surface trimmed in v only.
        need(args, 4, "trimv name surface v1 v2")?;
        // (Of an already trimmed surface, its basis keeping the other direction's trim, as
        // OCCT's.)
        let surface = self.surface_var(args[1])?;
        let (v1, v2) = (self.number(args[2])?, self.number(args[3])?);
        let (u1, u2, _, _) = kentron::geometry::ParametricSurface::bounds(&surface);
        let surface = match surface {
            Surface::Trimmed(t) => t.basis().clone(),
            s => s,
        };
        self.set(args[0], DrawObject::Surface(Surface::Trimmed(TrimmedSurface::new(surface, u1, u2, v1, v2))));
        Ok(String::new())
    }
    /// DRAW command `trinfo`: 8 uses in 6 converted scripts.
    ///
    /// Example: `trinfo b`
    pub fn trinfo(&mut self, args: &[&str]) -> DrawResult {
        // `trinfo shape`: the numbers of faces, of those without a mesh, and of triangles and
        // nodes of the faces' meshes, and their largest deflection, as OCCT reports them.
        need(args, 1, "trinfo shape")?;
        let shape = self.shape_var(args[0])?;
        let faces = self.model.unique_subshapes(shape, NodeKind::Face);
        let (mut triangles, mut nodes, mut empty, mut deflection) = (0usize, 0usize, 0usize, 0.0f64);
        for &face in &faces {
            match self.model.face_mesh(face) {
                Some(mesh) => {
                    triangles += mesh.triangles().len();
                    nodes += mesh.nodes().len();
                    deflection = deflection.max(mesh.deflection());
                }
                None => empty += 1,
            }
        }
        let pad = "                    ";
        let mut out = format!("\nThis shape contains {} faces.\n", faces.len());
        if empty > 0 {
            out += &format!("{pad}{empty} empty faces.\n");
        }
        out += &format!("{pad}{triangles} triangles.\n{pad}{nodes} nodes.\nMaximal deflection {}\n", format_real(deflection));
        Ok(out)
    }
    /// DRAW command `countshapes`: `nbshapes`'s report.
    ///
    /// Example: `countshapes result`
    pub fn countshapes(&mut self, args: &[&str]) -> DrawResult {
        self.nbshapes(args)
    }
    /// DRAW command `vsetdispmode`: a shape shown shaded (mode 1) is meshed, as OCCT's viewer
    /// meshes it (a deflection of 0.001 of its bounding box's largest side, 20 degrees),
    /// unless it has a mesh; otherwise nothing.
    ///
    /// Example: `vsetdispmode res 1`
    pub fn vsetdispmode(&mut self, args: &[&str]) -> DrawResult {
        let [name, mode] = args else { return Ok(String::new()) };
        let Ok(shape) = self.shape_var(name) else { return Ok(String::new()) };
        if self.number(mode).ok() != Some(1.0) {
            return Ok(String::new());
        }
        let faces = self.model.unique_subshapes(shape, NodeKind::Face);
        if faces.is_empty() || faces.iter().all(|&f| self.model.face_mesh(f).is_some()) {
            return Ok(String::new());
        }
        let b = kentron::properties::bounding_box(&self.model, shape);
        let size = b.max - b.min;
        let deflection = 0.001 * size.x.max(size.y).max(size.z);
        if deflection > 0.0 {
            kentron::meshing::mesh_with(&mut self.model, shape, kentron::meshing::MeshOptions { deflection, angle: 20f64.to_radians() })?;
        }
        Ok(String::new())
    }
    /// DRAW command `triarea`: the area of a shape's mesh, and of its faces.
    ///
    /// Example: `triarea res 0`
    pub fn triarea(&mut self, args: &[&str]) -> DrawResult {
        // `triarea shape [eps]`: the area of the triangles of the faces' meshes, and the faces'
        // area (`eps`, OCCT's precision for that, not taken into account).
        if args.is_empty() || args.len() > 2 {
            return Err(DrawError::Failed("usage: triarea shape [eps]".into()));
        }
        let shape = self.shape_var(args[0])?;
        let mut mesh_area = 0.0;
        for face in self.model.unique_subshapes(shape, NodeKind::Face) {
            if let Some(mesh) = self.model.face_mesh(face) {
                let p = mesh.nodes();
                for t in mesh.triangles() {
                    let (a, b, c) = (p[t[0] as usize], p[t[1] as usize], p[t[2] as usize]);
                    mesh_area += 0.5 * (b - a).cross(c - a).length();
                }
            }
        }
        let area = kentron::properties::area(&self.model, shape).mass;
        Ok(format!("{} {}", format_real(mesh_area), format_real(area)))
    }
    /// DRAW command `trotate`: 1438 uses in 1062 converted scripts.
    ///
    /// Example: `trotate s 0 0 0 0 0 1 90`
    pub fn trotate(&mut self, args: &[&str]) -> DrawResult {
        // `trotate name ... x y z dx dy dz angle`: moves the shapes (angle in degrees).
        // (OCCT's `-copy` and `-copymesh` change nothing here.)
        let args = &without_copy_flags(args);
        need(args, 8, "trotate name ... x y z dx dy dz angle")?;
        let n = args.len();
        let v = args[n - 7..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let direction = Direction::try_from_vector(Vector::new(v[3], v[4], v[5]))
            .ok_or_else(|| DrawError::Failed("trotate: null axis direction".into()))?;
        let t = Transform::rotation(Axis::new(Point::new(v[0], v[1], v[2]), direction), v[6].to_radians());
        self.move_shapes(&args[..n - 7], &t)
    }
    /// DRAW command `tscale`: 320 uses in 290 converted scripts.
    ///
    /// Example: `tscale s 0 0 0 SCALE1`
    pub fn tscale(&mut self, args: &[&str]) -> DrawResult {
        // `tscale name ... x y z factor`: scales the shapes about the point.
        // (OCCT's `-copy` and `-copymesh` change nothing here.)
        let options = transform_options(args);
        let args = &without_copy_flags(args);
        need(args, 5, "tscale name ... x y z factor")?;
        let n = args.len();
        let v = args[n - 4..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        // Like OCCT, a scaling rewrites the geometry (a placement keeps distances).
        let t = Transform::scaling(Point::new(v[0], v[1], v[2]), v[3]);
        for name in &args[..n - 4] {
            let shape = self.shape_var(name)?;
            let scaled = modify::transform(&mut self.model, shape, &t, options)?;
            self.set_shape_var(name, scaled)?;
        }
        Ok(String::new())
    }
    /// DRAW command `ttranslate`: 1310 uses in 910 converted scripts.
    ///
    /// Example: `ttranslate s2 0 7 7`
    pub fn ttranslate(&mut self, args: &[&str]) -> DrawResult {
        // `ttranslate name ... dx dy dz`: moves the shapes.
        // (OCCT's `-copy` and `-copymesh` change nothing here.)
        let args = &without_copy_flags(args);
        need(args, 4, "ttranslate name ... dx dy dz")?;
        let n = args.len();
        let v = args[n - 3..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let t = Transform::translation(Vector::new(v[0], v[1], v[2]));
        self.move_shapes(&args[..n - 3], &t)
    }
    /// DRAW command `tuyau`: 4 uses in 3 converted scripts.
    ///
    /// Example: `tuyau t curve c1 c2`
    pub fn tuyau(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("tuyau", args)
    }
    /// DRAW command `uiso`: 5 uses in 5 converted scripts.
    ///
    /// Example: `uiso l c pi`
    pub fn uiso(&mut self, args: &[&str]) -> DrawResult {
        // `uiso name surface u`: the iso-curve at u.
        need(args, 3, "uiso name surface u")?;
        let s = self.surface_var(args[1])?;
        let u = self.number(args[2])?;
        self.set(args[0], DrawObject::Curve(kentron::geometry::ParametricSurface::u_iso(&s, u)));
        Ok(String::new())
    }
    /// DRAW command `unifysamedom`: 144 uses in 124 converted scripts.
    ///
    /// Example: `unifysamedom result result +b`
    pub fn unifysamedom(&mut self, args: &[&str]) -> DrawResult {
        // `unifysamedom result shape [-f] [-e] [+b] [+i] [-nosafe] [-t tol] [-a angle]`: `-f`
        // keeps the faces, `-e` the edges; the other options (B-spline merging, internal
        // edges, tolerances) are not taken into account.
        need(args, 2, "unifysamedom result shape [options]")?;
        let shape = self.shape_var(args[1])?;
        let mut options = kentron::unify::UnifyOptions::default();
        let mut rest = args[2..].iter();
        while let Some(&a) = rest.next() {
            match a {
                "-f" => options.faces = false,
                "-e" => options.edges = false,
                "-t" | "-a" => {
                    rest.next();
                }
                _ => {}
            }
        }
        let result = kentron::unify::unify_same_domain_with(&mut self.model, shape, options)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `updatetolerance`: 29 uses in 13 converted scripts.
    ///
    /// Example: `updatetolerance a 1`
    pub fn updatetolerance(&mut self, args: &[&str]) -> DrawResult {
        // `updatetolerance [result] shape [verify_faces]`: makes the tolerances consistent, the
        // shape named `result`. (A last argument naming no shape is the flag, as OCCT reads it.)
        need(args, 1, "updatetolerance [result] shape [param]")?;
        let flag = args.len() >= 2 && self.shape_var(args[args.len() - 1]).is_err();
        let input = if flag { args[args.len() - 2] } else { args[args.len() - 1] };
        let shape = self.shape_var(input)?;
        modify::update_tolerances(&mut self.model, shape)?;
        if args[0] != input {
            self.set_shape_var(args[0], shape)?;
        }
        Ok(String::new())
    }
    /// DRAW command `updatevol`: 103 uses in 103 converted scripts.
    ///
    /// Example: `updatevol s_5 0 2 1 4 2 2`
    pub fn updatevol(&mut self, args: &[&str]) -> DrawResult {
        // `updatevol edge u1 r1 [u2 r2 ...]`: an edge of the evolving fillet, its radius through
        // the (u, r) pairs.
        need(args, 3, "updatevol edge u1 r1 [u2 r2 ...]")?;
        if (args.len() - 1) % 2 != 0 {
            return Err(DrawError::Failed("updatevol: parameters and radii come in pairs".into()));
        }
        let edge = self.shape_var(args[0])?;
        let values = self.numbers(&args[1..])?;
        let law: Vec<(f64, f64)> = values.chunks(2).map(|c| (c[0], c[1])).collect();
        let Some((_, _, edges)) = &mut self.evolving else {
            return Err(DrawError::Failed("updatevol: no mkevol before".into()));
        };
        edges.push((edge, law));
        Ok(String::new())
    }
    /// DRAW command `v2d`: 13 uses in 13 converted scripts.
    ///
    /// Example: `v2d`
    pub fn v2d(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("v2d", args)
    }
    /// DRAW command `vertex`: 1091 uses in 263 converted scripts.
    ///
    /// Example: `vertex va 0 0 0`
    pub fn vertex(&mut self, args: &[&str]) -> DrawResult {
        // `vertex name x y z`, or `vertex name point`; the forms on a curve or edge are not
        // supported.
        if args.len() == 2 {
            if let Some(DrawObject::Point(p)) = self.get(args[1]).cloned() {
                let vertex = make::vertex(&mut self.model, p);
                return self.set_shape_var(args[0], vertex);
            }
        }
        if args.len() != 4 {
            return self.not_implemented("vertex", args);
        }
        let p = Point::new(self.number(args[1])?, self.number(args[2])?, self.number(args[3])?);
        let vertex = make::vertex(&mut self.model, p);
        self.set_shape_var(args[0], vertex)
    }
    /// DRAW command `viso`: 3 uses in 3 converted scripts.
    ///
    /// Example: `viso l c 0`
    pub fn viso(&mut self, args: &[&str]) -> DrawResult {
        // `viso name surface v`: the iso-curve at v.
        need(args, 3, "viso name surface v")?;
        let s = self.surface_var(args[1])?;
        let v = self.number(args[2])?;
        self.set(args[0], DrawObject::Curve(kentron::geometry::ParametricSurface::v_iso(&s, v)));
        Ok(String::new())
    }
    /// DRAW command `wedge`: 13 uses in 13 converted scripts.
    ///
    /// Example: `wedge s 40 100 100 20 20 70 70`
    pub fn wedge(&mut self, args: &[&str]) -> DrawResult {
        // `wedge name [x y z nx ny nz vx vy vz] dx dy dz ltx | xmin zmin xmax zmax`: a box
        // whose top face shrinks (to a width ltx, or to the given rectangle).
        let values = args[1..].iter().map(|a| self.number(a)).collect::<DrawResult<Vec<f64>>>()?;
        let (frame, values) = match values.len() {
            4 | 7 => (Frame::WORLD, &values[..]),
            13 | 16 => {
                let v = &values;
                let normal = Direction::try_from_vector(Vector::new(v[3], v[4], v[5]))
                    .ok_or_else(|| DrawError::Failed("wedge: null direction".into()))?;
                let x = Direction::try_from_vector(Vector::new(v[6], v[7], v[8]))
                    .ok_or_else(|| DrawError::Failed("wedge: null direction".into()))?;
                (Frame::new(Point::new(v[0], v[1], v[2]), normal, x), &values[9..])
            }
            _ => return Err(DrawError::Failed("wedge: syntax error".into())),
        };
        let wedge = match *values {
            [dx, dy, dz, ltx] => Wedge::new(dx, dy, dz, ltx),
            [dx, dy, dz, x_min, z_min, x_max, z_max] => Wedge::with_top(dx, dy, dz, x_min, z_min, x_max, z_max),
            _ => unreachable!(),
        };
        let shape = wedge.with_frame(frame).build(&mut self.model)?;
        self.set_shape_var(args[0], shape)
    }
    /// DRAW command `wexplo`: 1 uses in 1 converted scripts.
    ///
    /// Example: `wexplo result_1 result`
    pub fn wexplo(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("wexplo", args)
    }
    /// DRAW command `wire`: 676 uses in 454 converted scripts.
    ///
    /// Example: `wire w w1 w2`
    pub fn wire(&mut self, args: &[&str]) -> DrawResult {
        // `wire name [-unsorted] e1|w1 e2|w2 ...`: the edges (and the edges of wires) in
        // order, or with -unsorted chained end to end whatever their order and orientation.
        // In order, as OCCT's BRepBuilderAPI_MakeWire takes them, each may also join the wire
        // so far at its start.
        need(args, 2, "wire name [-unsorted] edge|wire ...")?;
        let unsorted = args[1] == "-unsorted";
        let names = if unsorted { &args[2..] } else { &args[1..] };
        if names.iter().any(|a| a.starts_with('-')) {
            return self.not_implemented("wire", args);
        }
        let mut edges = Vec::new();
        for name in names {
            let shape = self.shape_var(name)?;
            match shape.kind() {
                NodeKind::Edge => edges.push(shape),
                NodeKind::Wire => edges.extend(self.model.explore(shape, NodeKind::Edge)),
                kind => return Err(DrawError::Failed(format!("{name} is a {kind:?}, not an edge or a wire"))),
            }
        }
        if unsorted {
            edges = self.chained(edges);
        }
        let wire = match make::wire(&mut self.model, &edges) {
            Ok(wire) => wire,
            Err(error) => {
                let chained = self.chained(edges.clone());
                if chained.len() != edges.len() {
                    return Err(error.into());
                }
                make::wire(&mut self.model, &chained)?
            }
        };
        self.set_shape_var(args[0], wire)
    }
    /// DRAW command `xdistcs`: 5 uses in 4 converted scripts.
    ///
    /// Example: `xdistcs ii s1 u1 u2 10 2.0e-7`
    pub fn xdistcs(&mut self, args: &[&str]) -> DrawResult {
        // `xdistcs curve surface t1 t2 count [tolerance [warning]]`: the distances to the
        // surface from points along the curve, an error where above the tolerance.
        need(args, 5, "xdistcs curve surface t1 t2 count [tolerance [warning]]")?;
        let curve = self.curve_var(args[0])?;
        let surface = self.surface_var(args[1])?;
        let (t1, t2) = (self.number(args[2])?, self.number(args[3])?);
        let count = (self.number(args[4])? as usize).max(2);
        let error = if args.len() > 5 { self.number(args[5])? } else { f64::INFINITY };
        let warning = if args.len() > 6 { self.number(args[6])? } else { f64::INFINITY };
        let mut out = String::new();
        let (mut max, mut at) = (0.0f64, 0.0);
        for i in 0..count {
            let t = if i + 1 == count { t2 } else { t1 + (t2 - t1) * i as f64 / (count - 1) as f64 };
            let d = kentron::geometry::extrema::nearest_on_surface(curve.point(t), &surface, None).2;
            if d > error {
                out += &format!("Error in {}:", args[0]);
            } else if d > warning {
                out += "Attention (critical value of tolerance) :";
            }
            out += &format!(" T={}\tD={}\n", format_real(t), format_real(d));
            if d > max {
                (max, at) = (d, t);
            }
        }
        out += &format!("Max distance = {}\nParam = {}\n", format_real(max), format_real(at));
        Ok(out)
    }
}

// Shims appended by tools/draw2rs.py. Implement in place; the tool never
// modifies an existing method.
#[allow(non_snake_case)]
impl Draw {
    /// DRAW command `DT_ClosedSplit`: 59 uses in 59 converted scripts.
    ///
    /// Example: `DT_ClosedSplit result a`
    pub fn DT_ClosedSplit(&mut self, args: &[&str]) -> DrawResult {
        // `DT_ClosedSplit result shape`: closed faces split in two. The result is a copy.
        let [result, shape] = args else { return self.not_implemented("DT_ClosedSplit", args) };
        let shape = self.shape_var(shape)?;
        let copy = kentron::modify::copy(&mut self.model, shape, false).unwrap_or(shape);
        let split = kentron::heal::split_closed_faces(&mut self.model, copy)?;
        self.set_shape_var(result, split)
    }
    /// DRAW command `DT_ToBspl`: 58 uses in 58 converted scripts.
    ///
    /// Example: `DT_ToBspl result a`
    pub fn DT_ToBspl(&mut self, args: &[&str]) -> DrawResult {
        // `DT_ToBspl result shape [options]`: the letters of the surfaces converted, `e`
        // (extrusions), `r` (revolutions and elementary surfaces), `o` (offsets), `p`
        // (planes); `ero` by default. The result is a copy.
        if !(2..=3).contains(&args.len()) {
            return self.not_implemented("DT_ToBspl", args);
        }
        let shape = self.shape_var(args[1])?;
        let options = args.get(2).copied().unwrap_or("ero");
        let kinds = kentron::heal::BSplineConversion {
            extrusion: options.contains('e'),
            revolution: options.contains('r'),
            offset: options.contains('o'),
            plane: options.contains('p'),
        };
        let copy = kentron::modify::copy(&mut self.model, shape, false).unwrap_or(shape);
        let result = kentron::heal::convert_to_bspline(&mut self.model, copy, &kinds)?;
        self.set_shape_var(args[0], result)
    }
    /// DRAW command `fixwgaps`: 12 uses in 12 converted scripts.
    ///
    /// Example: `fixwgaps result a 0.001`
    pub fn fixwgaps(&mut self, args: &[&str]) -> DrawResult {
        // `fixwgaps result shape [precision]`: the gaps in the wires covered (OCCT's
        // ShapeFix_Wireframe::FixWireGaps, which here grows tolerances as far as the gaps
        // need), on a copy.
        need(args, 2, "fixwgaps result shape [precision]")?;
        let shape = self.shape_var(args[1])?;
        let copy = modify::copy(&mut self.model, shape, false)?;
        let mut options = kentron::heal::FixOptions { max_tolerance: f64::INFINITY, ..Default::default() };
        if let Some(p) = args.get(2) {
            options.precision = self.number(p)?;
        }
        let result = kentron::heal::fix_shape(&mut self.model, copy, &options)?;
        self.set_shape_var(args[0], result)
    }
}

impl Draw {
    /// `edges` reordered and turned to run end to end (the chain grown from the first at
    /// both ends; those that do not join it follow as they were).
    fn chained(&self, edges: Vec<Shape>) -> Vec<Shape> {
        // (Start and end as the edge runs.)
        let ends = |e: Shape| {
            let (a, b) = self.model.edge_vertices(Shape { orientation: Orientation::Forward, ..e });
            let (a, b) = (self.model.point(a), self.model.point(b));
            if e.orientation == Orientation::Reversed { (b, a) } else { (a, b) }
        };
        let near = |p: Point, q: Point| p.distance(q) <= 1e-6;
        let mut rest = edges;
        if rest.is_empty() {
            return rest;
        }
        let mut chain = std::collections::VecDeque::from([rest.remove(0)]);
        loop {
            let (head, tail) = (ends(chain[0]).0, ends(chain[chain.len() - 1]).1);
            let Some(k) = rest.iter().position(|&e| {
                let (a, b) = ends(e);
                near(a, tail) || near(b, tail) || near(a, head) || near(b, head)
            }) else {
                break;
            };
            let e = rest.remove(k);
            let (a, b) = ends(e);
            if near(a, tail) {
                chain.push_back(e);
            } else if near(b, tail) {
                chain.push_back(e.reversed());
            } else if near(b, head) {
                chain.push_front(e);
            } else {
                chain.push_front(e.reversed());
            }
        }
        chain.into_iter().chain(rest).collect()
    }
}

// Shims appended by tools/draw2rs.py. Implement in place; the tool never
// modifies an existing method.
#[allow(non_snake_case)]
impl Draw {
    /// DRAW command `2dpscale`: 5 uses in 5 converted scripts.
    ///
    /// Example: `2dpscale result 0 0 0 1.1`
    pub fn _2dpscale(&mut self, args: &[&str]) -> DrawResult {
        // `2dpscale name ... x y factor`: the 2D curves and points scaled about (x, y) (names
        // that are not variables are passed over, as DRAW does).
        need(args, 4, "2dpscale name ... x y factor")?;
        let n = args.len();
        let v = self.numbers(&args[n - 3..])?;
        let names: Vec<&str> = args[..n - 3].iter().copied().filter(|a| self.get(a).is_some()).collect();
        self.transform_vars2(&names, &Transform2::scaling(Point2::new(v[0], v[1]), v[2]))
    }
    /// DRAW command `u`: 1 uses in 1 converted scripts.
    ///
    /// Example: `u`
    pub fn u(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("u", args)
    }
}

impl Draw {
    /// `fillet2d`'s and `chamfer2d`'s shared forms: on two edges (the corner edge as the
    /// result, the edges cut back) or on a wire of two edges (the wire rebuilt); `sizes`
    /// numbers after the shapes (a radius, or two lengths).
    fn corner_2d(&mut self, command: &str, args: &[&str], sizes: usize) -> DrawResult {
        let (shapes, numbers) = match args.len() {
            n if n == 3 + sizes => (&args[1..3], &args[3..]),
            n if n == 2 + sizes => (&args[1..2], &args[2..]),
            _ => return self.not_implemented(command, args),
        };
        let v = self.numbers(numbers)?;
        let edges: Vec<Shape> = match shapes {
            [e1, e2] => vec![self.shape_var(e1)?, self.shape_var(e2)?],
            [w] => {
                let w = self.shape_var(w)?;
                self.model.ordered_edges(w)
            }
            _ => unreachable!(),
        };
        let [e1, e2] = edges[..] else {
            return Err(DrawError::Failed(format!("{command}: the wire must have two edges")));
        };
        let corner = if sizes == 1 { blend::fillet_2d(&mut self.model, e1, e2, v[0])? } else { blend::chamfer_2d(&mut self.model, e1, e2, v[0], v[1])? };
        if let [n1, n2] = shapes {
            self.set_shape_var(n1, corner.first)?;
            self.set_shape_var(n2, corner.second)?;
            self.set_shape_var(args[0], corner.corner)
        } else {
            let wire = make::wire(&mut self.model, &[corner.first, corner.corner, corner.second])?;
            self.set_shape_var(args[0], wire)
        }
    }
}

/// `shape` and its sub-shapes numbered as OCCT's boolean data structure numbers them (its
/// `BOPDS_DS`, whose numbers `bopcheck` names shapes by): the shape, then each child and its
/// sub-shapes in turn, depth first, each once (a wire's children being its edges).
fn occt_numbering(model: &kentron::topology::Model, shape: Shape) -> Vec<Shape> {
    fn visit(model: &kentron::topology::Model, s: Shape, out: &mut Vec<Shape>, seen: &mut std::collections::HashSet<NumberingKey>) {
        let s = if s.kind() == NodeKind::CoEdge { model.explore(s, NodeKind::Edge).next().unwrap_or(s) } else { s };
        if !seen.insert(numbering_key(s)) {
            return;
        }
        out.push(s);
        for child in model.children(s).collect::<Vec<_>>() {
            visit(model, child, out, seen);
        }
    }
    let mut out = Vec::new();
    visit(model, shape, &mut out, &mut Default::default());
    out
}

/// A shape's node and placement (what [`occt_numbering`] tells shapes apart by), hashable.
type NumberingKey = (kentron::topology::NodeId, [[u64; 4]; 3]);

fn numbering_key(s: Shape) -> NumberingKey {
    // (`+ 0.0`: a negative zero as a zero, as `==` takes it.)
    (s.node, s.placement.rows().map(|row| row.map(|v| (v + 0.0).to_bits())))
}

/// DRAW's names of the kinds of shapes, in `nbshapes`' order.
const SHAPE_KINDS: [(&str, NodeKind); 8] = [
    ("vertex", NodeKind::Vertex),
    ("edge", NodeKind::Edge),
    ("wire", NodeKind::Wire),
    ("face", NodeKind::Face),
    ("shell", NodeKind::Shell),
    ("solid", NodeKind::Solid),
    ("compsolid", NodeKind::CompSolid),
    ("compound", NodeKind::Compound),
];

/// How many distinct sub-shapes of `kind` `shape` has (the shape itself included), as OCCT
/// counts them: each placement of a sub-shape apart, or with `shared` each TShape once.
fn count_subshapes(model: &kentron::topology::Model, shape: Shape, kind: NodeKind, shared: bool) -> usize {
    // (Compounds within compounds count too, as OCCT maps every sub-shape.)
    let found: Vec<Shape> = if kind == NodeKind::Compound {
        occt_numbering(model, shape).into_iter().filter(|s| s.kind() == kind).collect()
    } else {
        model.unique_subshapes(shape, kind)
    };
    if shared {
        let mut nodes: Vec<_> = found.iter().map(|s| s.node).collect();
        nodes.sort_unstable();
        nodes.dedup();
        nodes.len()
    } else {
        found.len()
    }
}

/// The report of `lprops`, `sprops` and `vprops`: the mass and the centre of gravity.
fn mass_report(p: properties::MassProperties) -> String {
    let c = p.center_of_mass;
    format!("Mass : {}\n\nCenter of gravity : \nX = {}\nY = {}\nZ = {}\n", format_real(p.mass), format_real(c.x), format_real(c.y), format_real(c.z))
}

/// `edges` in connected chains (through their shared vertices), each in order along it, the
/// edges oriented along it.
fn chain_edges(model: &kentron::topology::Model, edges: &[Shape]) -> Vec<Vec<Shape>> {
    let ends = |e: Shape| {
        let (a, b) = model.edge_vertices(e);
        (a.node, b.node)
    };
    let mut left: Vec<Shape> = edges.to_vec();
    let mut chains = Vec::new();
    while let Some(first) = left.pop() {
        let mut chain = std::collections::VecDeque::from([first]);
        loop {
            let (head, _) = ends(chain[0]);
            let (_, tail) = ends(chain[chain.len() - 1]);
            if let Some(k) = left.iter().position(|&e| { let (a, b) = ends(e); a == tail || b == tail }) {
                let e = left.remove(k);
                chain.push_back(if ends(e).0 == tail { e } else { e.reversed() });
            } else if let Some(k) = left.iter().position(|&e| { let (a, b) = ends(e); a == head || b == head }) {
                let e = left.remove(k);
                chain.push_front(if ends(e).1 == head { e } else { e.reversed() });
            } else {
                break;
            }
        }
        chains.push(chain.into_iter().collect());
    }
    chains
}

/// A chain of edges as one curve: an edge's curve over its range, or for several a curve
/// interpolated through points along them.
fn chain_curve(model: &kentron::topology::Model, chain: &[Shape]) -> DrawResult<Curve> {
    if let [e] = chain {
        let (c, first, last) = model.edge_curve(*e).ok_or_else(|| DrawError::Failed("an edge without a curve".into()))?;
        let (first, last) = if e.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
        return Ok(if first < last { Curve::Trimmed(TrimmedCurve::new(c, first, last)) } else { Curve::Trimmed(TrimmedCurve::new(c, last, first)) });
    }
    let mut samples = Vec::new();
    for &e in chain {
        let Some((c, first, last)) = model.edge_curve(e) else { continue };
        let (from, to) = if e.orientation == Orientation::Reversed { (last, first) } else { (first, last) };
        let n = 32;
        for i in 0..=n {
            if i == 0 && !samples.is_empty() {
                continue;
            }
            samples.push(c.point(from + (to - from) * i as f64 / n as f64));
        }
    }
    samples.dedup_by(|a, b| a.distance(*b) <= 1e-9);
    Ok(Curve::BSpline(kentron::geometry::fit::interpolate(&samples, false, 1e-12)?))
}

/// Where two surfaces meet in closed form (elementary ones, an extrusion along lines, or a
/// cylinder or cone and a quadric along branches made B-splines): the curves and isolated
/// points; `None` if not.
fn elementary_intersection(s1: &Surface, s2: &Surface, tolerance: f64) -> Option<(Vec<Curve>, Vec<Point>)> {
    use kentron::geometry::intersect::{extrusion_lines, quadrics, ruled_quadric, QuadricIntersection as Q};
    let lines = match (s1, s2) {
        (Surface::Extrusion(e), other) | (other, Surface::Extrusion(e)) => extrusion_lines(e, other, tolerance.max(1e-7)),
        _ => None,
    };
    if let Some(lines) = lines {
        return Some((lines.into_iter().map(Curve::Line).collect(), Vec::new()));
    }
    let elementary = |s: &Surface| matches!(s, Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_) | Surface::Sphere(_) | Surface::Torus(_));
    if !elementary(s1) || !elementary(s2) {
        return None;
    }
    let tolerance = tolerance.max(1e-7);
    Some(match quadrics(s1, s2, tolerance) {
        // (Coinciding surfaces are left to the general method, failing as OCCT's.)
        Q::Same => return None,
        Q::Empty => (Vec::new(), Vec::new()),
        Q::Points(p) => (Vec::new(), p),
        Q::Lines(c) => (c.into_iter().map(Curve::Line).collect(), Vec::new()),
        Q::Circles(c) => (c.into_iter().map(Curve::Circle).collect(), Vec::new()),
        Q::Ellipses(c) => (c.into_iter().map(Curve::Ellipse).collect(), Vec::new()),
        Q::Parabolas(c) => (c.into_iter().map(Curve::Parabola).collect(), Vec::new()),
        Q::Hyperbolas(c) => (c.into_iter().map(Curve::Hyperbola).collect(), Vec::new()),
        Q::NotAnalytic => {
            let branches = ruled_quadric(s1, s2, tolerance, None).or_else(|| ruled_quadric(s2, s1, tolerance, None))?;
            let mut curves = Vec::new();
            for branch in &branches.curves {
                let (a, b) = branch.domain();
                let n = 256;
                let points: Vec<Point> = (0..=n).map(|i| branch.point_clamped(a + (b - a) * i as f64 / n as f64)).collect();
                curves.push(fit_through_on_surfaces(&points, s1, s2)?);
            }
            (curves, Vec::new())
        }
    })
}

/// A curve through points along where two surfaces meet, moved exactly onto both, with
/// more of them where it bends. `None` for fewer than two distinct points.
fn fit_through_on_surfaces(points: &[Point], s1: &Surface, s2: &Surface) -> Option<Curve> {
    use kentron::geometry::intersect::surfaces_point;
    // (On the untrimmed surfaces: a curve's ends are on their trimming boundaries.)
    let untrimmed = |s: &Surface| match s {
        Surface::Trimmed(t) => t.basis().clone(),
        s => s.clone(),
    };
    let (s1, s2) = (&untrimmed(s1), &untrimmed(s2));
    let onto = |p: Point| surfaces_point(s1, s2, p, None).map(|(q, _)| q).filter(|q| q.distance(p) <= 1e-2 * (1.0 + (p - Point::ORIGIN).length()));
    let mut refined: Vec<Point> = points.iter().map(|&p| onto(p).unwrap_or(p)).collect();
    refined.dedup_by(|a, b| a.distance(*b) <= 1e-9);
    // Halve the steps where the curve turns by more than about a degree, or a step is
    // long for the curve.
    let length: f64 = refined.windows(2).map(|w| w[0].distance(w[1])).sum();
    for _ in 0..6 {
        let mut next = vec![refined[0]];
        let mut split = false;
        for i in 1..refined.len() {
            let (a, b) = (refined[i - 1], refined[i]);
            let turn = |p: Point, q: Point, r: Point| (q - p).angle(r - q);
            let bent = (i >= 2 && turn(refined[i - 2], a, b) > 0.02) || (i + 1 < refined.len() && turn(a, b, refined[i + 1]) > 0.02);
            if (bent || a.distance(b) > length / 128.0) && refined.len() < 4096 {
                if let Some(m) = onto(a + (b - a) * 0.5) {
                    if m.distance(a) > 1e-9 && m.distance(b) > 1e-9 {
                        next.push(m);
                        split = true;
                    }
                }
            }
            next.push(b);
        }
        refined = next;
        if !split {
            break;
        }
    }
    if refined.len() < 2 {
        return None;
    }
    kentron::geometry::fit::interpolate(&refined, false, 1e-12).ok().map(Curve::BSpline)
}

impl Draw {
    /// The helix commands: `helix`, `comphelix` and `spiral` (by heights, and pitches or
    /// numbers of turns as each PF is 1 or 0) and their `2` forms (by pitches and numbers of
    /// turns), with one diameter, one each end of each part, or the first and last.
    fn helix_command(&mut self, command: &str, args: &[&str]) -> DrawResult {
        use kentron::helix::{helix, HelixOptions, HelixPart};
        let usage = || DrawError::Failed(format!("{command}: wrong arguments"));
        need(args, 4, command)?;
        let n = self.number(args[1])? as usize;
        let v = self.numbers(&args[2..])?;
        let by_turns = command.ends_with('2');
        let diameters_given = match command.trim_end_matches('2') {
            "helix" => 1,
            "spiral" => 2,
            _ => n + 1,
        };
        let per_part = if by_turns { 2 } else { 3 };
        if v.len() != diameters_given + per_part * n || n == 0 {
            return Err(usage());
        }
        let (diameters, rest) = v.split_at(diameters_given);
        let column = |k: usize| &rest[k * n..(k + 1) * n];
        // Heights and pitches.
        let (heights, pitches): (Vec<f64>, Vec<f64>) = if by_turns {
            (column(0).iter().zip(column(1)).map(|(p, t)| p * t).collect(), column(0).to_vec())
        } else {
            let (h, p, is_pitch) = (column(0), column(1), column(2));
            (h.to_vec(), (0..n).map(|i| if is_pitch[i] != 0.0 { p[i] } else { h[i] / p[i] }).collect())
        };
        // The diameter at each part's end.
        let ends: Vec<f64> = match diameters_given {
            1 => vec![diameters[0]; n],
            2 => {
                let total: f64 = heights.iter().sum();
                let slope = (diameters[1] - diameters[0]) / total;
                let mut risen = 0.0;
                (0..n)
                    .map(|i| {
                        risen += heights[i];
                        if i + 1 == n { diameters[1] } else { diameters[0] + slope * risen }
                    })
                    .collect()
            }
            _ => diameters[1..].to_vec(),
        };
        let parts: Vec<HelixPart> = (0..n).map(|i| HelixPart { height: heights[i], pitch: pitches[i], end_diameter: ends[i] }).collect();
        let frame = self.helix_frame.unwrap_or(Frame::WORLD);
        let options = HelixOptions::default();
        match helix(&mut self.model, &frame, diameters[0], &parts, &options) {
            Ok(made) => {
                self.set_shape_var(args[0], made.wire)?;
                let warning = i32::from(made.tolerance_reached > options.tolerance);
                Ok(format!("WarningStatus = {warning}\nToleranceReached = {}\n", format_real(made.tolerance_reached)))
            }
            Err(e) => Ok(format!("ErrorStatus = 1 ({e})\n")),
        }
    }
}

// Shims appended by tools/draw2rs.py. Implement in place; the tool never
// modifies an existing method.
#[allow(non_snake_case)]
impl Draw {
    /// DRAW command `xdistef`: 2 uses in 1 converted scripts.
    ///
    /// Example: `xdistef e_0 f_0`
    pub fn xdistef(&mut self, args: &[&str]) -> DrawResult {
        self.not_implemented("xdistef", args)
    }
}
