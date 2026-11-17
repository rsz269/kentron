//! OCCT's arc length along a curve and its inverse: `CPnts_AbscissaPoint` (Gauss integration
//! of the speed `|C'(u)|`, and the parameter at a given arc length by `math_FunctionRoot` on
//! that integral) and `GCPnts_AbscissaPoint` (the same with the per-curve-type shortcuts:
//! lines, circles and two-pole polynomial Bezier and B-spline curves are parametrized
//! proportionally to length, and curves of several CN spans are taken span by span), with
//! the `math_GaussSingleIntegration` they integrate with. 3D curves only (OCCT's 2D
//! overloads are the same code over `Adaptor2d_Curve2d`).
//!
//! kentron's own arc length is [`measure`](super::measure); this port is for algorithms
//! ported from OCCT that must see OCCT's numbers (`Approx_CurvlinFunc` among them).
//!
//! The integrand, the speed of the curve, always evaluates, so `math_GaussSingleIntegration`
//! is always done and the `Standard_ConstructionError` OCCT raises when it is not cannot
//! happen; the port's lengths are plain numbers.

use super::adaptor::{CurveAdaptor, CurveKind, hunt};
use super::curve::Curve;
use super::function_root::function_root;
use super::function_roots::FunctionWithDerivative;
use super::{Continuity, ParametricCurve};
use crate::{Error, Result};

/// `Precision::Confusion()`.
const CONFUSION: f64 = 1e-7;
/// `Precision::PConfusion()`.
const P_CONFUSION: f64 = 1e-9;

/// `math::GaussPointsMax()`.
const GAUSS_POINTS_MAX: usize = 61;

// `Point[]` and `Weight[]` of OCCT's math.cxx, the rows for 1 to 24 points (the most
// `CPnts_AbscissaPoint` asks for), digits as OCCT writes them. Row `N` holds the
// `(N + 1) / 2` non-negative points in decreasing order (0 last for odd `N`) and their
// weights; index 0 is unused, as in OCCT.
#[allow(clippy::excessive_precision)]
#[rustfmt::skip]
const GAUSS_POINTS: [f64; 157] = [
    0.0e+00,
    // N = 1
    0.0e+00,
    // N = 2
    0.577350269189625764509148780501958e+00,
    // N = 3
    0.774596669241483377035853079956480e+00, 0.0e+00,
    // N = 4
    0.861136311594052575223946488892810e+00, 0.339981043584856264802665759103245e+00,
    // N = 5
    0.906179845938663992797626878299393e+00, 0.538469310105683091036314420700209e+00, 0.0e+00,
    // N = 6
    0.932469514203152027812301554493994e+00, 0.661209386466264513661399595019906e+00, 0.238619186083196908630501721680712e+00,
    // N = 7
    0.949107912342758524526189684047851e+00, 0.741531185599394439863864773280789e+00, 0.405845151377397166906606412076961e+00, 0.0e+00,
    // N = 8
    0.960289856497536231683560868569472e+00, 0.796666477413626739591553936475831e+00, 0.525532409916328985817739049189246e+00, 0.183434642495649804939476142360184e+00,
    // N = 9
    0.968160239507626089835576202903672e+00, 0.836031107326635794299429788069735e+00, 0.613371432700590397308702039341474e+00, 0.324253423403808929038538014643337e+00, 0.0e+00,
    // N = 10
    0.973906528517171720077964012084452e+00, 0.865063366688984510732096688423493e+00, 0.679409568299024406234327365114874e+00, 0.433395394129247190799265943165784e+00, 0.148874338981631210884826001129720e+00,
    // N = 11
    0.978228658146056992803938001122856e+00, 0.887062599768095299075157769303928e+00, 0.730152005574049324093416252031154e+00, 0.519096129206811815925725669458610e+00, 0.269543155952344972331531985400861e+00, 0.0e+00,
    // N = 12
    0.981560634246719250690549090149283e+00, 0.904117256370474856678465866119096e+00, 0.769902674194304687036893833212819e+00, 0.587317954286617447296702418940535e+00, 0.367831498998180193752691536643721e+00, 0.125233408511468915472441369463853e+00,
    // N = 13
    0.984183054718588149472829448807119e+00, 0.917598399222977965206547836500699e+00, 0.801578090733309912794206489582858e+00, 0.642349339440340220643984606995517e+00, 0.448492751036446852877912852127640e+00, 0.230458315955134794065528121097989e+00, 0.0e+00,
    // N = 14
    0.986283808696812338841597266704064e+00, 0.928434883663573517336391139377848e+00, 0.827201315069764993189794742650405e+00, 0.687292904811685470148019803019327e+00, 0.515248636358154091965290718551190e+00, 0.319112368927889760435671824168475e+00, 0.108054948707343662066244650219835e+00,
    // N = 15
    0.987992518020485428489565718586659e+00, 0.937273392400705904307758947710112e+00, 0.848206583410427216200648320774263e+00, 0.724417731360170047416186054613928e+00, 0.570972172608538847537226737253913e+00, 0.394151347077563369897207370981046e+00, 0.201194093997434522300628303394596e+00, 0.0e+00,
    // N = 16
    0.989400934991649932596154173450404e+00, 0.944575023073232576077988415534452e+00, 0.865631202387831743880467897712467e+00, 0.755404408355003033895101194847403e+00, 0.617876244402643748446671764048800e+00, 0.458016777657227386342419442983577e+00, 0.281603550779258913230460501460496e+00, 0.950125098376374401853193354249581e-01,
    // N = 17
    0.990575475314417335675434019940707e+00, 0.950675521768767761222716957895783e+00, 0.880239153726985902122955694488141e+00, 0.781514003896801406925230055520479e+00, 0.657671159216690765850302216643003e+00, 0.512690537086476967886246568629553e+00, 0.351231763453876315297185517095346e+00, 0.178484181495847855850677493654066e+00, 0.0e+00,
    // N = 18
    0.991565168420930946730016004706406e+00, 0.955823949571397755181195892928630e+00, 0.892602466497555739206060591127977e+00, 0.803704958972523115682417455014383e+00, 0.691687043060353207874891081288922e+00, 0.559770831073947534607871548525322e+00, 0.411751161462842646035931793833052e+00, 0.251886225691505509588972854877911e+00, 0.847750130417353012422618529357838e-01,
    // N = 19
    0.992406843843584403189017670254354e+00, 0.960208152134830030852778840685415e+00, 0.903155903614817901642660928533756e+00, 0.822714656537142824978922486712007e+00, 0.720966177335229378617095860824125e+00, 0.600545304661681023469638164946193e+00, 0.464570741375960945717267148104109e+00, 0.316564099963629831990117328849845e+00, 0.160358645640225375868096115740744e+00, 0.0e+00,
    // N = 20
    0.993128599185094924786122388473126e+00, 0.963971927277913791267666131192107e+00, 0.912234428251325905867752441208091e+00, 0.839116971822218823394529061698667e+00, 0.746331906460150792614305070356382e+00, 0.636053680726515025452836696226146e+00, 0.510867001950827098004364050955266e+00, 0.373706088715419560672548177024927e+00, 0.227785851141645078080496195368575e+00, 0.765265211334973337546404093988382e-01,
    // N = 21
    0.993752170620389500260242035940332e+00, 0.967226838566306294316622214904022e+00, 0.920099334150400828790187133718838e+00, 0.853363364583317283647250638585182e+00, 0.768439963475677908615877851306936e+00, 0.667138804197412319305966669990138e+00, 0.551618835887219807059018796724313e+00, 0.424342120207438783573668888543788e+00, 0.288021316802401096600792516064600e+00, 0.145561854160895090937030982338686e+00, 0.0e+00,
    // N = 22
    0.994294585482399292073031421168066e+00, 0.970060497835428727123950986748560e+00, 0.926956772187174000520692939267243e+00, 0.865812577720300136536425637014974e+00, 0.787816805979208162004277955410719e+00, 0.694487263186682780050689835761995e+00, 0.587640403506911592958876927638605e+00, 0.469355837986757026406330710966410e+00, 0.341935820892084225158147420427379e+00, 0.207860426688221285478846533919546e+00, 0.697392733197222212138417961186281e-01,
    // N = 23
    0.994769334997552123523925715452405e+00, 0.972542471218115231956024076799993e+00, 0.932971086826016102349196989085222e+00, 0.876752358270441667378156885913797e+00, 0.804888401618839892151118407002211e+00, 0.718661363131950194461624483747842e+00, 0.619609875763646156385097311649803e+00, 0.509501477846007549689793047866828e+00, 0.390301038030290831421488872880605e+00, 0.264135680970344930533869538283309e+00, 0.133256824298466110931742682241766e+00, 0.0e+00,
    // N = 24
    0.995187219997021360179997409695963e+00, 0.974728555971309498198391993040647e+00, 0.938274552002732758523649001703114e+00, 0.886415527004401034213154341985133e+00, 0.820001985973902921953949872662897e+00, 0.740124191578554364243828103101925e+00, 0.648093651936975569252495786909868e+00, 0.545421471388839535658375617218405e+00, 0.433793507626045138487084231913343e+00, 0.315042679696163374386793291319810e+00, 0.191118867473616309158639820757070e+00, 0.640568928626056260850430826247450e-01,
];
#[allow(clippy::excessive_precision)]
#[rustfmt::skip]
const GAUSS_WEIGHTS: [f64; 157] = [
    0.0e+00,
    // N = 1
    0.200000000000000000000000000000000e+01,
    // N = 2
    0.100000000000000000000000000000000e+01,
    // N = 3
    0.555555555555555555555555555555555e+00, 0.888888888888888888888888888888889e+00,
    // N = 4
    0.347854845137453857373063949221999e+00, 0.652145154862546142626936050778001e+00,
    // N = 5
    0.236926885056189087514264040719917e+00, 0.478628670499366468041291514835638e+00, 0.568888888888888888888888888888889e+00,
    // N = 6
    0.171324492379170345040296142172733e+00, 0.360761573048138607569833513837716e+00, 0.467913934572691047389870343989551e+00,
    // N = 7
    0.129484966168869693270611432679083e+00, 0.279705391489276667901467771423778e+00, 0.381830050505118944950369775488978e+00, 0.417959183673469387755102040816324e+00,
    // N = 8
    0.101228536290376259152531354309963e+00, 0.222381034453374470544355994426239e+00, 0.313706645877887287337962201986603e+00, 0.362683783378361982965150449277195e+00,
    // N = 9
    0.812743883615744119718921581105244e-01, 0.180648160694857404058472031242912e+00, 0.260610696402935462318742869418636e+00, 0.312347077040002840068630406584444e+00, 0.330239355001259763164525069286978e+00,
    // N = 10
    0.666713443086881375935688098933317e-01, 0.149451349150580593145776339657693e+00, 0.219086362515982043995534934228166e+00, 0.269266719309996355091226921569470e+00, 0.295524224714752870173892994651339e+00,
    // N = 11
    0.556685671161736664827537204425497e-01, 0.125580369464904624634694299223937e+00, 0.186290210927734251426097641431661e+00, 0.233193764591990479918523704843169e+00, 0.262804544510246662180688869890515e+00, 0.272925086777900630714483528336339e+00,
    // N = 12
    0.471753363865118271946159614850133e-01, 0.106939325995318430960254718194013e+00, 0.160078328543346226334652529543347e+00, 0.203167426723065921749064455809812e+00, 0.233492536538354808760849898924877e+00, 0.249147045813402785000562436042951e+00,
    // N = 13
    0.404840047653158795200215922009913e-01, 0.921214998377284479144217759538088e-01, 0.138873510219787238463601776868882e+00, 0.178145980761945738280046691996084e+00, 0.207816047536888502312523219306048e+00, 0.226283180262897238412090186039811e+00, 0.232551553230873910194589515268787e+00,
    // N = 14
    0.351194603317518630318328761381690e-01, 0.801580871597602098056332770628825e-01, 0.121518570687903184689414809072460e+00, 0.157203167158193534569601938623896e+00, 0.185538397477937813741716590125112e+00, 0.205198463721295603965924065661222e+00, 0.215263853463157790195876443316267e+00,
    // N = 15
    0.307532419961172683546283935771943e-01, 0.703660474881081247092674164506892e-01, 0.107159220467171935011869546685734e+00, 0.139570677926154314447804794511208e+00, 0.166269205816993933553200860480891e+00, 0.186161000015562211026800561866565e+00, 0.198431485327111576456118326443407e+00, 0.202578241925561272880620199967656e+00,
    // N = 16
    0.271524594117540948517805724560368e-01, 0.622535239386478928628438369945628e-01, 0.951585116824927848099251076018181e-01, 0.124628971255533872052476282192097e+00, 0.149595988816576732081501730547117e+00, 0.169156519395002538189312079030817e+00, 0.182603415044923588866763667969201e+00, 0.189450610455068496285396723208291e+00,
    // N = 17
    0.241483028685479319601100262874223e-01, 0.554595293739872011294401653581783e-01, 0.850361483171791808835353701908650e-01, 0.111883847193403971094788385627117e+00, 0.135136368468525473286319981701965e+00, 0.154045761076810288081431594801964e+00, 0.168004102156450044509970663788286e+00, 0.176562705366992646325270990114256e+00, 0.179446470356206525458265644261812e+00,
    // N = 18
    0.216160135264833103133427102655130e-01, 0.497145488949697964533349462010895e-01, 0.764257302548890565291296776163320e-01, 0.100942044106287165562813984928046e+00, 0.122555206711478460184519126798015e+00, 0.140642914670650651204731303753862e+00, 0.154684675126265244925418003835190e+00, 0.164276483745832722986053776466578e+00, 0.169142382963143591840656470134736e+00,
    // N = 19
    0.194617882297264770363120414644793e-01, 0.448142267656996003328381574039431e-01, 0.690445427376412265807082579987094e-01, 0.914900216224499994644620941319334e-01, 0.111566645547333994716023901674946e+00, 0.128753962539336227675515784862361e+00, 0.142606702173606611775746109431414e+00, 0.152766042065859666778855400902497e+00, 0.158968843393954347649956439460513e+00, 0.161054449848783695979163625324809e+00,
    // N = 20
    0.176140071391521183118619623491302e-01, 0.406014298003869413310399522706322e-01, 0.626720483341090635695065351819669e-01, 0.832767415767047487247581432329156e-01, 0.101930119817240435036750135466348e+00, 0.118194531961518417312377377718604e+00, 0.131688638449176626898494499748159e+00, 0.142096109318382051329298325072027e+00, 0.149172986472603746787828736999336e+00, 0.152753387130725850698084331956079e+00,
    // N = 21
    0.160172282577743333242246168579398e-01, 0.369537897708524937999506683027812e-01, 0.571344254268572082836358264526385e-01, 0.761001136283793020170516533126966e-01, 0.934444234560338615532897410907539e-01, 0.108797299167148377663474578081143e+00, 0.121831416053728534195367177105957e+00, 0.132268938633337461781052574507564e+00, 0.139887394791073154722133423860591e+00, 0.144524403989970059063827166570068e+00, 0.146081133649690427191985147672526e+00,
    // N = 22
    0.146279952982722006849910980409279e-01, 0.337749015848141547933022468735675e-01, 0.522933351526832859403120512609446e-01, 0.697964684245204880949614189253340e-01, 0.859416062170677274144436813628141e-01, 0.100414144442880964932078837856422e+00, 0.112932296080539218393400607417128e+00, 0.123252376810512424285560986159639e+00, 0.131173504787062370732964992516422e+00, 0.136541498346015171352573831229592e+00, 0.139251872855631993375410248339356e+00,
    // N = 23
    0.134118594871417720813094934295282e-01, 0.309880058569794443106942196749048e-01, 0.480376717310846685716410715037491e-01, 0.642324214085258521271696153694344e-01, 0.792814117767189549228925246444880e-01, 0.929157660600351474770186174421938e-01, 0.104892091464541410074086184918183e+00, 0.114996640222411364941643512970691e+00, 0.123049084306729530467578400635032e+00, 0.128905722188082149978595339377180e+00, 0.132462039404696617371642464587308e+00, 0.133654572186106175351457110560180e+00,
    // N = 24
    0.123412297999871995468056670241292e-01, 0.285313886289336631813078161424755e-01, 0.442774388174198061686027483473371e-01, 0.592985849154367807463677585232428e-01, 0.733464814110803057340336151500117e-01, 0.861901615319532759171852029475615e-01, 0.976186521041138882698806644404639e-01, 0.107444270115965634782577342387576e+00, 0.115505668053725601353344483876487e+00, 0.121670472927803391204463153443133e+00, 0.125837456346828296121375382504667e+00, 0.127938195346752156974056165213338e+00,
];

/// `math_GaussSingleIntegration::Perform(F, Lower, Upper, Order)`: Gauss-Legendre quadrature
/// of `f` over `[lower, upper]` with `order` points, summed by symmetric pairs.
fn gauss_perform(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize) -> f64 {
    assert!((1..=24).contains(&order), "Gauss points beyond 24 are not tabulated here (CPnts_AbscissaPoint never asks for them)");
    // math::GaussPoints and math::GaussWeights: the row of `order` points starts after
    // those of 1 to order - 1.
    let start: usize = (1..order).map(|i| (i + 1) >> 1).sum();
    let point = |j: usize| GAUSS_POINTS[start + j];
    let weight = |j: usize| GAUSS_WEIGHTS[start + j];
    // The change of variable to [lower, upper].
    let xm = 0.5 * (upper + lower);
    let xr = 0.5 * (upper - lower);
    let mut val = 0.0;
    let ind = order / 2;
    let ind1 = order.div_ceil(2);
    if ind1 > ind {
        // An odd count: the middle point.
        val = f(xm);
        val *= weight(ind1);
    }
    // The sum over the Gauss points, by symmetry.
    for j in 1..=ind {
        let dx = xr * point(j);
        let f1 = f(xm - dx);
        let f2 = f(xm + dx);
        let ft = f1 + f2;
        val += weight(j) * ft;
    }
    // Back to [lower, upper].
    val *= xr;
    val
}

/// `math_GaussSingleIntegration(F, Lower, Upper, Order)`: the integral of `f` over
/// `[lower, upper]` (negative if `upper < lower`).
#[doc(alias = "math_GaussSingleIntegration")]
pub(crate) fn gauss_single_integration(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize) -> f64 {
    gauss_perform(f, lower, upper, GAUSS_POINTS_MAX.min(order))
}

/// `math_GaussSingleIntegration(F, Lower, Upper, Order, Tol)`: the integral over 1, 2, 4, ...
/// equal subintervals until two successive sums differ by at most `tolerance` (13 rounds at
/// most).
#[doc(alias = "math_GaussSingleIntegration")]
pub(crate) fn gauss_single_integration_with_tolerance(f: &dyn Fn(f64) -> f64, lower: f64, upper: f64, order: usize, tolerance: f64) -> f64 {
    let order = GAUSS_POINTS_MAX.min(order);
    // IterMax.
    const ITERATION_LIMIT: usize = 13;
    let mut iteration = 1;
    let mut interval_count: usize = 1;
    let mut len = gauss_perform(f, lower, upper, order);
    loop {
        let old_len = len;
        len = 0.0;
        interval_count *= 2;
        let du = (upper - lower) / interval_count as f64;
        for i in 1..=interval_count {
            len += gauss_perform(f, lower + (i - 1) as f64 * du, lower + i as f64 * du, order);
        }
        iteration += 1;
        if !((old_len - len).abs() > tolerance && iteration <= ITERATION_LIMIT) {
            break;
        }
    }
    len
}

/// `f3d`: the speed `|C'(x)|`.
fn speed(curve: &dyn CurveAdaptor, x: f64) -> f64 {
    curve.d1(x).1.length()
}

/// The curve behind an adaptor of a kind `curve_type` names (OCCT's `Bezier()`,
/// `BSpline()`, `Circle()`, which raise `Standard_NoSuchObject` otherwise).
fn basis(curve: &dyn CurveAdaptor) -> &Curve {
    curve.basis_curve().expect("the adaptor names its curve kind but gives no curve (Standard_NoSuchObject)")
}

/// `order(C)`: the Gauss point count for a kind of curve.
fn integration_order(curve: &dyn CurveAdaptor) -> usize {
    match curve.curve_type() {
        Some(CurveKind::Line) => 2,
        Some(CurveKind::Parabola) => 5,
        Some(CurveKind::Bezier) => match basis(curve) {
            Curve::Bezier(b) => 24.min(2 * b.degree()),
            _ => unreachable!("a Bezier adaptor over another curve"),
        },
        Some(CurveKind::BSpline) => match basis(curve) {
            Curve::BSpline(b) => 24.min(2 * b.poles().len() - 1),
            _ => unreachable!("a B-spline adaptor over another curve"),
        },
        _ => 10,
    }
}

/// `CPnts_AbscissaPoint::Length(C, U1, U2)` and `Length(C, U1, U2, Tol)`: the length of the
/// curve between `u1` and `u2` by Gauss integration of its speed, with `tolerance` by
/// repeated halving until the sums settle.
#[doc(alias = "CPnts_AbscissaPoint::Length")]
pub(crate) fn integrated_length(curve: &dyn CurveAdaptor, u1: f64, u2: f64, tolerance: Option<f64>) -> f64 {
    let f = |x: f64| speed(curve, x);
    let value = match tolerance {
        Some(tolerance) => gauss_single_integration_with_tolerance(&f, u1, u2, integration_order(curve), tolerance),
        None => gauss_single_integration(&f, u1, u2, integration_order(curve)),
    };
    value.abs()
}

/// `CPnts_MyRootFunction`: the arc length from `x0` minus the length sought, whose
/// derivative is the speed.
struct RootFunction<'a> {
    curve: &'a dyn CurveAdaptor,
    order: usize,
    x0: f64,
    length: f64,
    /// `myTol`: at most 0 for the plain integration.
    tolerance: f64,
}

impl RootFunction<'_> {
    fn integral(&self, x: f64) -> f64 {
        let f = |t: f64| speed(self.curve, t);
        if self.tolerance <= 0.0 {
            gauss_single_integration(&f, self.x0, x, self.order)
        } else {
            gauss_single_integration_with_tolerance(&f, self.x0, x, self.order, self.tolerance)
        }
    }
}

impl FunctionWithDerivative for RootFunction<'_> {
    fn value(&mut self, x: f64, value: &mut f64) -> bool {
        *value = self.integral(x) - self.length;
        true
    }

    fn derivative(&mut self, x: f64, derivative: &mut f64) -> bool {
        *derivative = speed(self.curve, x);
        true
    }

    fn values(&mut self, x: f64, value: &mut f64, derivative: &mut f64) -> bool {
        *value = self.integral(x) - self.length;
        *derivative = speed(self.curve, x);
        true
    }

    fn get_state_number(&mut self) -> i32 {
        0
    }
}

/// `CPnts_AbscissaPoint`: the parameter at a given arc length from a given parameter, by
/// `math_FunctionRoot` on the integrated length, searched within the initialized range
/// widened by its own length on either side.
#[doc(alias = "CPnts_AbscissaPoint")]
#[derive(Clone, Copy, Default)]
pub(crate) struct CurveAbscissaPoint<'a> {
    /// The curve and Gauss order of `myF`.
    curve: Option<&'a dyn CurveAdaptor>,
    order: usize,
    done: bool,
    /// `myL`: the length of the initialized range.
    length: f64,
    param: f64,
    u_min: f64,
    u_max: f64,
}

impl<'a> CurveAbscissaPoint<'a> {
    /// `CPnts_AbscissaPoint(C, Abscissa, U0, Resolution)`: `Init(C, Resolution)` then
    /// `Perform(Abscissa, U0, Resolution)`.
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, resolution: f64) -> Self {
        let mut point = Self::default();
        point.init_with_tolerance(curve, resolution);
        point.perform(abscissa, u0, resolution);
        point
    }

    /// `CPnts_AbscissaPoint(C, Abscissa, U0, Ui, Resolution)`: `Init(C)` then
    /// `Perform(Abscissa, U0, Ui, Resolution)`.
    pub(crate) fn with_guess(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, ui: f64, resolution: f64) -> Self {
        let mut point = Self::default();
        point.init(curve);
        point.perform_with_guess(abscissa, u0, ui, resolution);
        point
    }

    /// `Init(C)`: over the whole curve.
    pub(crate) fn init(&mut self, curve: &'a dyn CurveAdaptor) {
        self.init_range(curve, curve.first_parameter(), curve.last_parameter(), None);
    }

    /// `Init(C, Tol)`: over the whole curve, its length to `tolerance`.
    pub(crate) fn init_with_tolerance(&mut self, curve: &'a dyn CurveAdaptor, tolerance: f64) {
        self.init_range(curve, curve.first_parameter(), curve.last_parameter(), Some(tolerance));
    }

    /// `Init(C, U1, U2)` and `Init(C, U1, U2, Tol)`: the length between `u1` and `u2`, and
    /// the search range, `[u1, u2]` widened by its own length on either side.
    pub(crate) fn init_range(&mut self, curve: &'a dyn CurveAdaptor, u1: f64, u2: f64, tolerance: Option<f64>) {
        self.curve = Some(curve);
        self.order = integration_order(curve);
        self.length = integrated_length(curve, u1, u2, tolerance);
        self.u_min = u1.min(u2);
        self.u_max = u1.max(u2);
        let du = self.u_max - self.u_min;
        self.u_min -= du;
        self.u_max += du;
    }

    /// `Perform(Abscissa, U0, Resolution)`: guesses the parameter as if the speed were
    /// uniform over a third of the search range ("exercise: why 3?", asks OCCT).
    pub(crate) fn perform(&mut self, abscissa: f64, u0: f64, resolution: f64) {
        if self.length < CONFUSION {
            // Leave less violently: the parameter is not expected to move much.
            self.done = true;
            self.param = u0;
        } else {
            let ui = u0 + (abscissa / self.length) * (self.u_max - self.u_min) / 3.0;
            self.perform_with_guess(abscissa, u0, ui, resolution);
        }
    }

    /// `Perform(Abscissa, U0, Ui, Resolution)`: the root from the guess `ui`, the lengths
    /// integrated plainly.
    pub(crate) fn perform_with_guess(&mut self, abscissa: f64, u0: f64, ui: f64, resolution: f64) {
        self.solve(abscissa, u0, ui, resolution, -1.0);
    }

    /// `AdvPerform(Abscissa, U0, Ui, Resolution)`: the same, the lengths integrated to a
    /// tenth of `resolution`.
    pub(crate) fn adv_perform(&mut self, abscissa: f64, u0: f64, ui: f64, resolution: f64) {
        self.solve(abscissa, u0, ui, resolution, resolution / 10.0);
    }

    fn solve(&mut self, abscissa: f64, u0: f64, ui: f64, resolution: f64, integration_tolerance: f64) {
        if self.length < CONFUSION {
            // Leave less violently.
            self.done = true;
            self.param = u0;
            return;
        }
        self.done = false;
        let curve = self.curve.expect("CPnts_AbscissaPoint performed before Init");
        let mut function = RootFunction { curve, order: self.order, x0: u0, length: abscissa, tolerance: integration_tolerance };
        let solution = function_root(&mut function, ui, resolution, self.u_min, self.u_max, 100);
        // (OCCT suspends its check of the solution's accuracy, "temporarily" since 1997.)
        if solution.done {
            self.done = true;
            self.param = solution.root;
        }
    }

    pub(crate) fn is_done(&self) -> bool {
        self.done
    }

    /// `Parameter()`: `None` where OCCT raises `StdFail_NotDone`.
    pub(crate) fn parameter(&self) -> Option<f64> {
        self.done.then_some(self.param)
    }

    /// `SetParameter(P)`.
    pub(crate) fn set_parameter(&mut self, p: f64) {
        self.param = p;
        self.done = true;
    }
}

/// `GCPnts_AbscissaType`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum AbscissaType {
    /// `GCPnts_LengthParametrized`: length is the parameter times a ratio.
    LengthParametrized,
    /// `GCPnts_Parametrized`: one smooth span.
    Parametrized,
    /// `GCPnts_AbsComposite`: several CN spans.
    Composite,
}

/// `computeType`: the type, and for a length-parametrized curve the length per unit of
/// parameter.
fn compute_type(curve: &dyn CurveAdaptor) -> (AbscissaType, f64) {
    if curve.intervals(Continuity::CN).len() - 1 > 1 {
        return (AbscissaType::Composite, 1.0);
    }
    match curve.curve_type() {
        Some(CurveKind::Line) => (AbscissaType::LengthParametrized, 1.0),
        Some(CurveKind::Circle) => match basis(curve) {
            Curve::Circle(c) => (AbscissaType::LengthParametrized, c.radius()),
            _ => unreachable!("a circle adaptor over another curve"),
        },
        Some(CurveKind::Bezier) => match basis(curve) {
            Curve::Bezier(b) if b.poles().len() == 2 && !b.is_rational() => (AbscissaType::LengthParametrized, b.dn(0.0, 1).length()),
            _ => (AbscissaType::Parametrized, 1.0),
        },
        Some(CurveKind::BSpline) => match basis(curve) {
            Curve::BSpline(b) if b.poles().len() == 2 && !b.is_rational() => (AbscissaType::LengthParametrized, b.dn(b.first_parameter(), 1).length()),
            _ => (AbscissaType::Parametrized, 1.0),
        },
        _ => (AbscissaType::Parametrized, 1.0),
    }
}

/// `GCPnts_AbscissaPoint::Length(C, U1, U2)` and `Length(C, U1, U2, Tol)` (the `length`
/// template): exact for length-parametrized curves, else integrated, span by span for a
/// curve of several CN spans.
#[doc(alias = "GCPnts_AbscissaPoint::Length")]
pub(crate) fn length_between(curve: &dyn CurveAdaptor, u1: f64, u2: f64, tolerance: Option<f64>) -> f64 {
    let (kind, ratio) = compute_type(curve);
    match kind {
        AbscissaType::LengthParametrized => (u2 - u1).abs() * ratio,
        AbscissaType::Parametrized => integrated_length(curve, u1, u2, tolerance),
        AbscissaType::Composite => {
            let ti = curve.intervals(Continuity::CN);
            let uu1 = u1.min(u2);
            let uu2 = u1.max(u2);
            let mut l = 0.0;
            for k in 0..ti.len() - 1 {
                if ti[k] > uu2 {
                    break;
                }
                if ti[k + 1] < uu1 {
                    continue;
                }
                l += integrated_length(curve, ti[k].max(uu1), ti[k + 1].min(uu2), tolerance);
            }
            l
        }
    }
}

/// `GCPnts_AbscissaPoint::Length(C)` and `Length(C, Tol)`: the length of the whole curve.
#[doc(alias = "GCPnts_AbscissaPoint::Length")]
pub(crate) fn length(curve: &dyn CurveAdaptor, tolerance: Option<f64>) -> f64 {
    length_between(curve, curve.first_parameter(), curve.last_parameter(), tolerance)
}

/// `Compute`: the parameter at arc length `abscis` from `u0`, `ui` a guess, `epsilon` the
/// parameter resolution.
fn compute<'a>(computer: &mut CurveAbscissaPoint<'a>, curve: &'a dyn CurveAdaptor, mut abscis: f64, mut u0: f64, mut ui: f64, epsilon: f64) {
    // The easy solution.
    if abscis.abs() <= CONFUSION {
        computer.set_parameter(u0);
        return;
    }
    let (kind, ratio) = compute_type(curve);
    match kind {
        AbscissaType::LengthParametrized => computer.set_parameter(u0 + abscis / ratio),
        AbscissaType::Parametrized => {
            computer.init(curve);
            computer.perform_with_guess(abscis, u0, ui, epsilon);
        }
        AbscissaType::Composite => {
            let ti = curve.intervals(Continuity::CN);
            let interval_count = ti.len() - 1;
            // `aTI(i)`, 1-based.
            let t = |i: usize| ti[i - 1];
            let mut sign = 1.0;
            let mut index = hunt(&ti, u0);
            let mut direction = 1;
            if abscis < 0.0 {
                direction = 0;
                abscis = -abscis;
                sign = -1.0;
            }
            while index >= 1 && index <= interval_count {
                let l = integrated_length(curve, u0, t(index + direction), None);
                if (l - abscis).abs() <= CONFUSION {
                    computer.set_parameter(t(index + direction));
                    return;
                }
                if l > abscis {
                    if ui < t(index) || ui > t(index + 1) {
                        ui = (abscis / l) * (t(index + 1) - u0);
                        ui = if direction == 1 { u0 + ui } else { u0 - ui };
                    }
                    computer.init_range(curve, t(index), t(index + 1), None);
                    computer.perform_with_guess(sign * abscis, u0, ui, epsilon);
                    return;
                }
                u0 = t(index + direction);
                abscis -= l;
                if direction == 1 {
                    index += 1;
                } else {
                    index -= 1;
                }
            }
            // Push a little outside the limits ("hairy!!!", says OCCT).
            ui = u0 + 0.1;
            computer.init_range(curve, u0, u0 + 0.2, None);
            computer.perform_with_guess(sign * abscis, u0, ui, epsilon);
        }
    }
}

/// `AdvCompute`: [`compute`] with the lengths integrated to `epsilon` (used as both the
/// length tolerance and the parameter resolution).
fn adv_compute<'a>(computer: &mut CurveAbscissaPoint<'a>, curve: &'a dyn CurveAdaptor, mut abscis: f64, mut u0: f64, mut ui: f64, epsilon: f64) {
    let (kind, ratio) = compute_type(curve);
    match kind {
        AbscissaType::LengthParametrized => computer.set_parameter(u0 + abscis / ratio),
        AbscissaType::Parametrized => {
            computer.init_with_tolerance(curve, epsilon);
            computer.adv_perform(abscis, u0, ui, epsilon);
        }
        AbscissaType::Composite => {
            let ti = curve.intervals(Continuity::CN);
            let interval_count = ti.len() - 1;
            // `aTI(i)`, 1-based.
            let t = |i: usize| ti[i - 1];
            let mut sign = 1.0;
            let mut index = hunt(&ti, u0);
            let mut direction = 1;
            if abscis < 0.0 {
                direction = 0;
                abscis = -abscis;
                sign = -1.0;
            }

            // Starting before the first interval, forwards.
            if index == 0 && direction > 0 {
                let l = integrated_length(curve, u0, t(index + direction), Some(epsilon));
                if (l - abscis).abs() <= epsilon {
                    computer.set_parameter(t(index + direction));
                    return;
                }
                if l > abscis {
                    if ui > t(index + 1) {
                        ui = (abscis / l) * (t(index + 1) - u0);
                        ui += u0;
                    }
                    computer.init_range(curve, u0, t(index + 1), Some(epsilon));
                    computer.adv_perform(sign * abscis, u0, ui, epsilon);
                    return;
                }
                u0 = t(index + direction);
                abscis -= l;
                index += 1;
            }

            while index >= 1 && index <= interval_count {
                let l = integrated_length(curve, u0, t(index + direction), Some(epsilon));
                if (l - abscis).abs() <= P_CONFUSION {
                    computer.set_parameter(t(index + direction));
                    return;
                }
                if l > abscis {
                    if ui < t(index) || ui > t(index + 1) {
                        ui = (abscis / l) * (t(index + 1) - u0);
                        ui = if direction == 1 { u0 + ui } else { u0 - ui };
                    }
                    computer.init_range(curve, t(index), t(index + 1), Some(epsilon));
                    computer.adv_perform(sign * abscis, u0, ui, epsilon);
                    return;
                }
                u0 = t(index + direction);
                abscis -= l;
                if direction == 1 {
                    index += 1;
                } else {
                    index -= 1;
                }
            }

            // Push a little outside the limits ("hairy!!!"), but not past the ends of a
            // curve that is not periodic.
            let is_non_periodic = !curve.is_periodic();
            ui = u0 + sign * 0.1;
            let mut u1 = u0 + sign * 0.2;
            if is_non_periodic {
                if sign > 0.0 {
                    ui = ui.min(curve.last_parameter());
                    u1 = u1.min(curve.last_parameter());
                } else {
                    ui = ui.max(curve.first_parameter());
                    u1 = u1.max(curve.first_parameter());
                }
            }
            computer.init_range(curve, u0, u1, Some(epsilon));
            computer.adv_perform(sign * abscis, u0, ui, epsilon);
        }
    }
}

/// `GCPnts_AbscissaPoint`: the parameter of the point at a given arc length (signed: negative
/// goes backwards) from a given parameter.
#[doc(alias = "GCPnts_AbscissaPoint")]
#[derive(Clone, Copy)]
pub(crate) struct AbscissaPoint<'a> {
    computer: CurveAbscissaPoint<'a>,
}

impl<'a> AbscissaPoint<'a> {
    /// `GCPnts_AbscissaPoint(C, Abscissa, U0)`: the guess from the curve's average speed, the
    /// parameter resolution that of `Precision::Confusion()`. An error where OCCT raises
    /// `Standard_ConstructionError`: a curve shorter than `Precision::Confusion()`.
    pub(crate) fn new(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64) -> Result<Self> {
        let l = length(curve, None);
        if l < CONFUSION {
            return Err(Error::InvalidInput("GCPnts_AbscissaPoint: the curve is shorter than Precision::Confusion()".into()));
        }
        let ui = u0 + (abscissa / l) * (curve.last_parameter() - curve.first_parameter());
        let mut computer = CurveAbscissaPoint::default();
        compute(&mut computer, curve, abscissa, u0, ui, curve.resolution(CONFUSION));
        Ok(AbscissaPoint { computer })
    }

    /// `GCPnts_AbscissaPoint(Tol, C, Abscissa, U0)`: lengths to `tolerance`, which is also
    /// the parameter resolution.
    pub(crate) fn with_tolerance(tolerance: f64, curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64) -> Self {
        let l = length(curve, Some(tolerance));
        // (OCCT's check of a too short curve is commented out here.)
        let ui = if l >= CONFUSION { u0 + (abscissa / l) * (curve.last_parameter() - curve.first_parameter()) } else { u0 };
        let mut computer = CurveAbscissaPoint::default();
        adv_compute(&mut computer, curve, abscissa, u0, ui, tolerance);
        AbscissaPoint { computer }
    }

    /// `GCPnts_AbscissaPoint(C, Abscissa, U0, Ui)`: from the guess `ui`.
    pub(crate) fn with_guess(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, ui: f64) -> Self {
        let mut computer = CurveAbscissaPoint::default();
        compute(&mut computer, curve, abscissa, u0, ui, curve.resolution(CONFUSION));
        AbscissaPoint { computer }
    }

    /// `GCPnts_AbscissaPoint(C, Abscissa, U0, Ui, Tol)`: from the guess `ui`, lengths to
    /// `tolerance`, which is also the parameter resolution.
    pub(crate) fn with_guess_and_tolerance(curve: &'a dyn CurveAdaptor, abscissa: f64, u0: f64, ui: f64, tolerance: f64) -> Self {
        let mut computer = CurveAbscissaPoint::default();
        adv_compute(&mut computer, curve, abscissa, u0, ui, tolerance);
        AbscissaPoint { computer }
    }

    pub(crate) fn is_done(&self) -> bool {
        self.computer.is_done()
    }

    /// `Parameter()`: `None` where OCCT raises `StdFail_NotDone`.
    pub(crate) fn parameter(&self) -> Option<f64> {
        self.computer.parameter()
    }
}

#[cfg(test)]
mod tests;
