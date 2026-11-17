//! Tests against OCCT 8.1: the same curves and planes built in OCCT, run through
//! `IntCurveSurface_HInter::Perform(curve, GeomAdaptor_Surface(Geom_Plane(P, N)))`, the values
//! printed with `%.17g`.

#![allow(clippy::excessive_precision, clippy::approx_constant)]

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::bezier::BezierCurve;
use crate::geometry::bspline::BSplineCurve;
use crate::geometry::curve::{Circle, Ellipse, Hyperbola, Line, Parabola};
use crate::geometry::ParametricCurve;

/// The plane OCCT's `Geom_Plane(gp_Pnt(p), gp_Dir(n))` makes.
pub(crate) fn plane(p: [f64; 3], n: [f64; 3]) -> Plane {
    plane_through(Point::new(p[0], p[1], p[2]), Direction::new(n[0], n[1], n[2]))
}

fn pt(p: [f64; 3]) -> Point {
    Point::new(p[0], p[1], p[2])
}

/// `gp_Ax2(gp_Pnt(0.1, 0.2, 0.3), gp_Dir(0.2, 0.3, 1), gp_Dir(1, 0, 0))`, as OCCT 8.1
/// prints it.
fn frame_a1() -> Frame {
    Frame::from_axes_unchecked(
        Point::new(0.1, 0.2, 0.3),
        Direction::new(0.9821414205253256, -0.0540628304876326, -0.18020943495877534),
        Direction::new(0.0, 0.95782628522115154, -0.28734788556634544),
        Direction::new(0.18814417367671948, 0.2822162605150792, 0.94072086838359736),
    )
}

/// A cubic B-spline waving above and below z = 0 along X.
pub(crate) fn wavy_bspline() -> BSplineCurve {
    BSplineCurve::new(
        [[0.0, 0.0, 0.0], [1.0, 2.0, -1.0], [2.0, -1.0, 1.5], [3.0, 1.5, -2.0], [4.0, -2.0, 1.0], [5.0, 0.5, -0.5], [6.0, 1.0, 0.8]].into_iter().map(pt).collect(),
        vec![0.0, 0.3, 0.5, 0.8, 1.0],
        vec![4, 1, 1, 1, 4],
        3,
    )
}

/// A cubic B-spline through 17 poles turning twice round the Z axis, rising 0.2 a pole.
pub(crate) fn helix_bspline() -> BSplineCurve {
    let poles = (0..=16).map(|k| {
        let a = k as f64 * PI / 4.0;
        Point::new(2.0 * a.cos(), 2.0 * a.sin(), 0.2 * k as f64)
    });
    let mut mults = vec![1; 15];
    mults[0] = 4;
    mults[14] = 4;
    BSplineCurve::new(poles.collect(), (0..15).map(|k| k as f64).collect(), mults, 3)
}

fn periodic_bspline() -> BSplineCurve {
    let poles = (0..8).map(|k| {
        let t = k as f64 * PI / 4.0;
        Point::new(3.0 * t.cos(), 2.0 * t.sin(), 0.5 * (3.0 * t).sin())
    });
    BSplineCurve::new_periodic(poles.collect(), (0..9).map(|k| k as f64).collect(), vec![1; 9], 3)
}

/// A parabolic arc touching z = 0 at its middle.
pub(crate) fn tangent_bspline() -> BSplineCurve {
    BSplineCurve::new(vec![pt([-1.0, 0.0, 1.0]), pt([0.0, 0.0, -1.0]), pt([1.0, 0.0, 1.0])], vec![0.0, 1.0], vec![3, 3], 2)
}

pub(crate) fn whole(curve: Curve) -> GeomAdaptor {
    let (first, last) = (curve.first_parameter(), curve.last_parameter());
    GeomAdaptor::new(&curve, first, last)
}

/// An expected point: `w`, `u`, `v`, the point and OCCT's transition number
/// (0 tangent, 1 in, 2 out).
type Expected = (f64, f64, f64, [f64; 3], i32);

/// Compares with OCCT 8.1's points (same count, order and transitions); returns the
/// largest deviation of a parameter or a coordinate.
fn check(name: &str, curve: &dyn CurveAdaptor, plane: &Plane, parallel: bool, expected: &[Expected]) -> f64 {
    let got = intersect_curve_plane(curve, plane).unwrap();
    assert_eq!(got.parallel, parallel, "{name}: IsParallel");
    assert_eq!(got.points.len(), expected.len(), "{name}: NbPoints, got {:?}", got.points);
    let mut deviation: f64 = 0.0;
    for (g, e) in got.points.iter().zip(expected) {
        let transition = match g.transition {
            Transition::Tangent => 0,
            Transition::In => 1,
            Transition::Out => 2,
        };
        assert_eq!(transition, e.4, "{name}: transition at w = {}", g.w);
        for (a, b) in [(g.w, e.0), (g.u, e.1), (g.v, e.2), (g.point.x, e.3[0]), (g.point.y, e.3[1]), (g.point.z, e.3[2])] {
            deviation = deviation.max((a - b).abs());
        }
    }
    assert!(deviation < 1e-12, "{name}: deviation {deviation:e} from OCCT, got {:?}", got.points);
    deviation
}

#[test]
fn lines() {
    let l1 = Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::new(1.0, 1.0, 1.0)));
    let pl = plane([0.5, 0.2, 0.1], [0.3, -0.4, 1.0]);
    let mut d = check("line", &GeomAdaptor::new(&l1, -10.0, 10.0), &pl, false, &[(-4.4840870907060939, -0.61692117905487309, 2.1684097517365446, [-1.5888888888888895, -0.58888888888888946, 0.41111111111111054], 2)]);
    d = d.max(check("line_short", &GeomAdaptor::new(&l1, -1.0, 1.0), &pl, false, &[]));
    let l2 = Curve::Line(Line::new(Point::new(1.0, 2.0, 3.0), Direction::new(1.0, 1.0, 0.0)));
    d = d.max(check("line_parallel", &GeomAdaptor::new(&l2, -10.0, 10.0), &plane([0.0; 3], [0.0, 0.0, 1.0]), true, &[]));
    eprintln!("lines: max deviation from OCCT {d:e}");
}

#[test]
fn circles() {
    let c1 = Curve::Circle(Circle::new(frame_a1(), 2.0));
    let pl = plane([0.3, 0.1, 0.2], [1.0, 0.2, 0.1]);
    let first = (1.6519655924154593, 2.0493794546202562, -0.44570830990140126, [-0.059264376938268221, 2.1181122905853114, -0.24358081178793983], 1);
    let second = (4.9695660474371062, -1.804342837063839, 0.56682765390803258, [0.59961824201538005, -1.6801522186769005, 0.76412201719999406], 2);
    let mut d = check("circle", &whole(c1.clone()), &pl, false, &[first, second]);
    d = d.max(check("circle_trimmed", &GeomAdaptor::new(&c1, 1.0, 4.0), &pl, false, &[first]));
    d = d.max(check("circle_far", &whole(c1), &plane([5.0, 0.1, 0.2], [1.0, 0.2, 0.1]), false, &[]));
    let c0 = Curve::Circle(Circle::new(Frame::WORLD, 1.5));
    d = d.max(check("circle_in_plane", &whole(c0.clone()), &plane([0.0; 3], [0.0, 0.0, 1.0]), true, &[]));
    d = d.max(check("circle_parallel", &whole(c0.clone()), &plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]), true, &[]));
    d = d.max(check("circle_tangent", &whole(c0), &plane([1.5, 0.0, 0.0], [1.0, 0.0, 0.0]), false, &[(0.0, 0.0, 0.0, [1.5, 0.0, 0.0], 0)]));
    eprintln!("circles: max deviation from OCCT {d:e}");
}

#[test]
fn ellipses_parabolas_hyperbolas() {
    let e1 = Curve::Ellipse(Ellipse::new(frame_a1(), 3.0, 1.5));
    let pl = plane([0.3, 0.1, 0.2], [1.0, 0.2, 0.1]);
    let first = (1.596731242074831, 1.5647552178405828, -0.31837701004827129, [0.023593302512801562, 1.6404621384556122, -0.11685730203924399], 1);
    let second = (4.8568463527319432, -1.3629613367130484, 0.45085805126444523, [0.52415390470794654, -1.2451225094928788, 0.64870597190627421], 2);
    let mut d = check("ellipse", &whole(e1.clone()), &pl, false, &[first, second]);
    d = d.max(check("ellipse_trimmed", &GeomAdaptor::new(&e1, 0.5, 2.5), &pl, false, &[first]));
    d = d.max(check(
        "ellipse_tilted",
        &whole(e1),
        &plane([0.3, -0.4, 0.5], [0.3, 0.7, 0.2]),
        false,
        &[
            (2.8036967300048619, -3.2231287231607424, -0.17284165453279599, [-2.6798165227784101, 0.82930046237872923, 0.66717316584206321], 1),
            (5.372719809113752, 1.7271417128004392, 0.19767610584026674, [1.9072723232074598, -1.0342043813030311, 0.3088068497494173], 2),
        ],
    ));
    let pb = Curve::Parabola(Parabola::new(Frame::WORLD, 0.7));
    d = d.max(check(
        "parabola",
        &GeomAdaptor::new(&pb, -5.0, 5.0),
        &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        false,
        &[
            (1.6733200530681511, 0.0, -1.6733200530681511, [1.0000000000000002, 1.6733200530681511, 0.0], 2),
            (-1.6733200530681509, 0.0, 1.6733200530681509, [0.99999999999999989, -1.6733200530681509, 0.0], 1),
        ],
    ));
    d = d.max(check(
        "parabola_tilted",
        &GeomAdaptor::new(&pb, -5.0, 5.0),
        &plane([0.5, 0.3, 0.0], [1.0, 0.4, 0.2]),
        false,
        &[
            (-1.9916424134538624, -2.4681744150937037, -2.7755575615628914e-17, [1.4166569653815451, -1.9916424134538624, 0.0], 1),
            (0.87164241345386229, 0.61567772143943433, 2.7755575615628914e-17, [0.27134303461845494, 0.87164241345386229, 0.0], 2),
        ],
    ));
    let hy = Curve::Hyperbola(Hyperbola::new(Frame::WORLD, 2.0, 1.0));
    d = d.max(check(
        "hyperbola",
        &GeomAdaptor::new(&hy, -2.0, 2.0),
        &plane([3.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        false,
        &[
            (0.96242365011920694, 0.0, -1.1180339887498949, [3.0, 1.1180339887498949, 0.0], 2),
            (-0.96242365011920694, 0.0, 1.1180339887498949, [3.0, -1.1180339887498949, 0.0], 1),
        ],
    ));
    d = d.max(check(
        "hyperbola_tilted",
        &GeomAdaptor::new(&hy, -2.0, 2.0),
        &plane([3.0, 0.5, 0.0], [1.0, -0.3, 0.1]),
        false,
        &[
            (1.0591042946819693, 0.80234827597729774, 0.0, [3.2305530805898788, 1.2685102686329301, 0.0], 2),
            (-0.75682342280903592, -1.389782401798352, 0.0, [2.6006489654459264, -0.83117011518024508, 0.0], 1),
        ],
    ));
    eprintln!("ellipses, parabolas, hyperbolas: max deviation from OCCT {d:e}");
}

#[test]
fn bsplines() {
    let bs1 = Curve::BSpline(wavy_bspline());
    let a = whole(bs1.clone());
    let mut d = check(
        "bspline_wavy",
        &a,
        &plane([2.0, 0.0, 0.05], [0.1, 0.05, 1.0]),
        false,
        &[
            (0.19514912648254937, -0.48288083475384957, 0.5779050522902488, [1.5166582290579542, 0.57719114830889806, 0.069474619678759675], 2),
            (0.3640122033128817, 0.34898314351589477, 0.12092598030428257, [2.3466533026106386, 0.12077659670148165, 0.0092958399038621797], 1),
            (0.71688897022776665, 1.8323197274876555, -0.60500244207912446, [3.8262176349078989, -0.6042550638542491, -0.10240901029807736], 2),
        ],
    );
    d = d.max(check(
        "bspline_once",
        &a,
        &plane([2.5, 0.0, 0.0], [1.0, 0.0, 0.0]),
        false,
        &[(0.40141929031152657, -0.26809546901620362, -0.26596656928637741, [2.5, 0.26596656928637741, -0.26809546901620362], 2)],
    ));
    d = d.max(check(
        "bspline_twice",
        &a,
        &plane([0.0, 0.3, 0.0], [0.0, 1.0, 0.2]),
        false,
        &[
            (0.018753387333396338, -0.16240306275302613, 0.18290941185043155, [0.18290941185043155, 0.33184986100172742, -0.15924930500863727], 2),
            (0.24800835399535767, 0.28526969869368712, 1.809695378972404, [1.809695378972404, 0.24405400922015999, 0.27972995389920097], 1),
            (0.43548408759590962, -0.53592452451488071, 2.635492277529802, [2.635492277529802, 0.40510344647362739, -0.52551723236813674], 2),
            (0.53136790868262829, -0.925467555181678, 3.0147721292142036, [3.0147721292142036, 0.48149912011801466, -0.90749560059007384], 1),
            (0.92529447680502019, 0.0093077527404667598, 5.1016713763219528, [5.1016713763219528, 0.2981745995057179, 0.0091270024714109194], 2),
        ],
    ));
    d = d.max(check("bspline_none", &a, &plane([0.0, 0.0, 10.0], [0.0, 0.0, 1.0]), false, &[]));
    // Evaluated past its knots, as `GeomFill_GuideTrihedronPlan` trims its guide.
    d = d.max(check(
        "bspline_extended",
        &GeomAdaptor::new(&bs1, -0.2, 1.2),
        &plane([-0.1, 0.0, 0.0], [1.0, 0.0, 0.0]),
        false,
        &[(-0.0098691886864338794, 0.10695771281055244, 0.20991498585330834, [-0.10000000000000001, -0.20991498585330834, 0.10695771281055244], 2)],
    ));
    eprintln!("bsplines: max deviation from OCCT {d:e}");
}

#[test]
fn helix_and_periodic_bsplines() {
    let helix = whole(Curve::BSpline(helix_bspline()));
    let mut d = check(
        "helix_axis_plane",
        &helix,
        &plane([0.0; 3], [1.0, 0.1, 0.0]),
        false,
        &[
            (1.1722155177597275, 1.7687453790859866, 0.4249894302055584, [-0.17599674322026221, 1.759967432202622, 0.4249894302055584], 1),
            (5.1271535617471304, -1.8043390081901827, 1.2254307123494264, [0.17953844168958319, -1.7953844168958379, 1.2254307123494264], 2),
            (9.1271535617471304, 1.8043390081901831, 2.0254307123494262, [-0.1795384416895828, 1.7953844168958384, 2.0254307123494262], 1),
            (13.032150262477503, -1.7591583718706159, 2.8247586835696854, [0.17504280034805139, -1.7504280034805166, 2.8247586835696854], 2),
        ],
    );
    d = d.max(check(
        "helix_level",
        &helix,
        &plane([0.0, 0.0, 1.3], [0.0, 0.0, 1.0]),
        false,
        &[(5.4999999999999991, 0.68984788275433362, -1.6654401145198812, [0.68984788275433362, -1.6654401145198812, 1.3], 2)],
    ));
    let periodic = whole(Curve::BSpline(periodic_bspline()));
    d = d.max(check(
        "periodic_wave",
        &periodic,
        &plane([0.0, 0.0, 0.1], [0.0, 0.0, 1.0]),
        false,
        &[
            (0.15730109387909313, 1.6636942994544726, 1.4229717212638364, [1.6636942994544726, 1.4229717212638364, 0.10000000000000001], 1),
            (1.8426989061209067, -1.6636942994544717, 1.4229717212638364, [-1.6636942994544717, 1.4229717212638364, 0.099999999999999978], 2),
            (2.7021161652909123, -2.6314048994926296, 0.41762177082328089, [-2.6314048994926296, 0.41762177082328089, 0.10000000000000006], 1),
            (4.5780450654558713, -0.87954219805162626, -1.7047375480442886, [-0.87954219805162626, -1.7047375480442886, 0.1000000000000002], 2),
            (5.4219549345441296, 0.87954219805162692, -1.7047375480442881, [0.87954219805162692, -1.7047375480442881, 0.099999999999999867], 1),
            (7.2978838347090882, 2.6314048994926291, 0.41762177082328095, [2.6314048994926291, 0.41762177082328095, 0.1000000000000001], 2),
        ],
    ));
    d = d.max(check(
        "periodic_x",
        &periodic,
        &plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
        false,
        &[
            (0.5176380902050417, -0.071394991575998934, -1.6748598491509408, [0.99999999999999989, 1.6748598491509408, -0.071394991575998934], 1),
            (5.4823619097949585, 0.071394991575998684, 1.6748598491509405, [0.99999999999999944, -1.6748598491509405, 0.071394991575998684], 2),
        ],
    ));
    eprintln!("helix and periodic: max deviation from OCCT {d:e}");
}

#[test]
fn tangent_and_lying_in_the_plane() {
    let z0 = plane([0.0; 3], [0.0, 0.0, 1.0]);
    let mut d = check("tangent", &whole(Curve::BSpline(tangent_bspline())), &z0, false, &[(0.5, 0.0, 0.0, [0.0, 0.0, 0.0], 0)]);
    let flat = BSplineCurve::new(vec![pt([0.0, 0.0, 0.0]), pt([1.0, 2.0, 0.0]), pt([2.0, -1.0, 0.0]), pt([3.0, 0.0, 0.0])], vec![0.0, 1.0], vec![4, 4], 3);
    // A curve in the plane: one null interval, no points.
    d = d.max(check("in_plane", &whole(Curve::BSpline(flat)), &z0, false, &[]));
    let partly = BSplineCurve::new(
        vec![pt([0.0, 0.0, 1.0]), pt([1.0, 0.0, 0.0]), pt([2.0, 0.0, 0.0]), pt([3.0, 0.0, 0.0]), pt([4.0, 0.0, 0.0]), pt([5.0, 0.0, -1.0])],
        vec![0.0, 1.0, 2.0, 3.0],
        vec![3, 1, 2, 3],
        2,
    );
    d = d.max(check(
        "in_plane_partly",
        &whole(Curve::BSpline(partly)),
        &z0,
        false,
        &[(1.0909090909090908, 1.5950413223140494, 0.0, [1.5950413223140494, 0.0, 0.0], 0), (2.0, 3.0, 0.0, [3.0, 0.0, 0.0], 0)],
    ));
    let bezier = BezierCurve::new(vec![pt([0.0, 0.0, -1.0]), pt([1.0, 1.0, 2.0]), pt([2.0, 0.0, -2.0]), pt([3.0, 1.0, 2.0]), pt([4.0, 0.0, -1.0])]);
    d = d.max(check(
        "bezier",
        &whole(Curve::Bezier(bezier)),
        &plane([0.0, 0.0, 0.2], [0.0, 0.1, 1.0]),
        false,
        &[
            (0.20301556153520747, 0.43994981058533456, -0.81206224614082989, [0.81206224614082989, 0.43776642335824822, 0.15622335766417517], 2),
            (0.40406856243989003, 0.50181286357468691, -1.6162742497595601, [1.6162742497595601, 0.49932246178258499, 0.15006775382174209], 1),
            (0.5959314375601088, 0.50181286357468702, -2.3837257502404352, [2.3837257502404352, 0.4993224617825851, 0.15006775382174142], 2),
            (0.79698443846479283, 0.43994981058533433, -3.1879377538591713, [3.1879377538591713, 0.43776642335824817, 0.15622335766417694], 1),
        ],
    ));
    eprintln!("tangent and in-plane: max deviation from OCCT {d:e}");
}

/// `math_DirectPolynomialRoots` and `math_TrigonometricFunctionRoots`.
#[test]
fn closed_form_solvers() {
    use crate::geometry::polynomial_roots::{PolynomialRoots, trigonometric_roots};
    let cases: [(PolynomialRoots, &[f64]); 10] = [
        (PolynomialRoots::quartic(1.0, -10.0, 35.0, -50.0, 24.0), &[2.0000000000000062, 1.0, 4.000000000000008, 2.999999999999988]),
        (PolynomialRoots::quartic(2.0, 0.5, -7.0, 1.25, 0.3), &[-2.0698814860101558, -0.13598560060645851, 1.6286602601584022, 0.32720682645821192]),
        (PolynomialRoots::quartic(1.0, 0.0, -2.0, 0.0, 1.0), &[-1.0, -1.0, 1.0, 1.0]),
        (PolynomialRoots::quartic(0.3, 1.0, 2.0, 3.0, 4.0), &[]),
        (PolynomialRoots::cubic(1.0, -6.0, 11.0, -6.0), &[3.0, 2.0, 1.0]),
        (PolynomialRoots::cubic(1.0, 0.0, 0.0, -8.0), &[2.0]),
        (PolynomialRoots::cubic(2.0, -3.0, 0.0, 1.0), &[1.0, 1.0, -0.5]),
        (PolynomialRoots::quadratic(1.0, -3.0, 2.0), &[2.0, 1.0]),
        (PolynomialRoots::quadratic(1.0, 2.0, 1.0), &[-1.0, -1.0]),
        (PolynomialRoots::quadratic(1.0, 0.0, 1.0), &[]),
    ];
    let mut d: f64 = 0.0;
    for (i, (got, expected)) in cases.iter().enumerate() {
        assert!(got.done && !got.infinite, "polynomial {i}");
        assert_eq!(got.roots.len(), expected.len(), "polynomial {i}: {:?}", got.roots);
        for (a, b) in got.roots.iter().zip(expected.iter()) {
            d = d.max((a - b).abs());
        }
    }
    let two_pi = 2.0 * PI;
    let trig: [(_, &[f64]); 6] = [
        (trigonometric_roots(0.0, 0.0, 1.2, -0.7, 0.3, 0.0, two_pi), &[1.260381158954275, 3.9666552513725919]),
        (trigonometric_roots(1.0, 0.3, 0.2, -0.5, -0.4, 0.0, two_pi), &[0.90709972904880676, 2.838586163728817]),
        (trigonometric_roots(0.0, 0.0, 0.0, 2.0, 1.0, 0.0, two_pi), &[5.7595865315812871, 3.6651914291880923]),
        (trigonometric_roots(0.0, 0.0, 2.0, 0.0, -1.0, 0.0, two_pi), &[1.0471975511965979, 5.2359877559829879]),
        (trigonometric_roots(0.0, 0.5, 0.0, 0.3, 0.0, 0.0, two_pi), &[0.0, 3.1415926535897931, 1.8754889808102941, 4.4076963263692921]),
        (trigonometric_roots(2.0, -1.0, 0.5, 0.25, -1.0, -1.0, 3.0), &[0.59406222780050044, 1.946042804525558]),
    ];
    for (i, (got, expected)) in trig.iter().enumerate() {
        assert!(got.done && !got.infinite, "trigonometric {i}");
        assert_eq!(got.roots.len(), expected.len(), "trigonometric {i}: {:?}", got.roots);
        for (a, b) in got.roots.iter().zip(expected.iter()) {
            d = d.max((a - b).abs());
        }
    }
    assert!(d < 1e-13, "deviation {d:e} from OCCT");
    eprintln!("closed-form solvers: max deviation from OCCT {d:e}");
}
