//! Tests against OCCT 8.1: the same curves built in OCCT, run through the same calls
//! (`GCPnts_AbscissaPoint`, `CPnts_AbscissaPoint`, `Approx_CurvlinFunc`), the values printed
//! with `%.17g`.

use super::*;
use crate::geometry::adaptor::GeomAdaptor;
use crate::geometry::curve::{BSplineCurve, BezierCurve, Circle, Ellipse, Line};
use crate::geometry::curvlin::CurvlinFunc;
use crate::math::{Direction, Frame, Point};

/// The values, in the order the tests compute them (all 17 digits OCCT printed).
#[allow(clippy::excessive_precision)]
const OCCT_VALUES: &[(&str, f64)] = &[
    ("cubic_c1 length", 11.450122872778373),
    ("cubic_c1 length_tol", 11.450126309234733),
    ("cubic_c1 length_ab", 6.644327695496882),
    ("cubic_c1 length_ab_tol", 6.6443311402506549),
    ("cubic_c1 cpnts_ab", 6.6368613440288371),
    ("cubic_c1 cpnts_ab_tol", 6.6443311402506513),
    ("cubic_c1 abscissa_fwd", 1.6241599758988357),
    ("cubic_c1 abscissa_back", 3.453495002639603),
    ("cubic_c1 abscissa_mid", 3.2791751455328142),
    ("cubic_c1 abscissa_tol_fwd", 1.6241599834461886),
    ("cubic_c1 abscissa_tol_back", 3.1545024731314908),
    ("cubic_c1 abscissa_guess_tol", 2.5238523441049892),
    ("cubic_c1 abscissa_guess", 1.9681424260365978),
    ("cubic_c1 cpnts_abscissa", 1.3145721252980938),
    ("cubic_c1 curvlin_length", 11.450126309244368),
    ("cubic_c1 s_of_u0", 0.13312612629609363),
    ("cubic_c1 s_of_u1", 0.39901262472867854),
    ("cubic_c1 s_of_u2", 0.83047511594121615),
    ("cubic_c1 u_of_s0", 0.0),
    ("cubic_c1 u_of_s1", 0.099163487654757299),
    ("cubic_c1 u_of_s2", 1.6241606331410698),
    ("cubic_c1 u_of_s3", 2.0821362198758937),
    ("cubic_c1 u_of_s4", 2.1102700859044385),
    ("cubic_c1 u_of_s5", 3.8763209581070432),
    ("cubic_c1 u_of_s6", 4.0000000000000018),
    ("cubic_c1 u_of_s7", 0.77089653127532376),
    ("cubic_c1 nb_intervals_c2", 2.0),
    ("cubic_c1 interval_c2_1", 0.0),
    ("cubic_c1 interval_c2_2", 0.4646753793345878),
    ("cubic_c1 interval_c2_3", 1.0),
    ("cubic_c1 trim_first", 0.14999999999999999),
    ("cubic_c1 trim_last", 0.84999999999999998),
    ("cubic_c1 eval_s0_d0_0", 1.6503321597299114),
    ("cubic_c1 eval_s0_d0_1", 0.79060380478802506),
    ("cubic_c1 eval_s0_d0_2", 0.43329267529363069),
    ("cubic_c1 eval_s0_d1_0", 9.7178093840458892),
    ("cubic_c1 eval_s0_d1_1", -5.1022162795162913),
    ("cubic_c1 eval_s0_d1_2", 3.2614356210170183),
    ("cubic_c1 eval_s0_d2_0", 54.650700123545732),
    ("cubic_c1 eval_s0_d2_1", 82.033904912766872),
    ("cubic_c1 eval_s0_d2_2", -34.503321378445371),
    ("cubic_c1 eval_s1_d0_0", 3.6232043298138406),
    ("cubic_c1 eval_s1_d0_1", 1.1697569214622765),
    ("cubic_c1 eval_s1_d0_2", 1.197750391392538),
    ("cubic_c1 eval_s1_d1_0", 4.0184184359574111),
    ("cubic_c1 eval_s1_d1_1", -9.5232324735435299),
    ("cubic_c1 eval_s1_d1_2", 4.9260277126761638),
    ("cubic_c1 eval_s1_d2_0", 37.457735024333616),
    ("cubic_c1 eval_s1_d2_1", -9.7235814252099146),
    ("cubic_c1 eval_s1_d2_2", -49.354326357740618),
    ("cubic_c1 eval_s2_d0_0", 5.7268742640344881),
    ("cubic_c1 eval_s2_d0_1", 1.4470806505935776),
    ("cubic_c1 eval_s2_d0_2", -0.10990203478745803),
    ("cubic_c1 eval_s2_d1_0", 10.220366426911671),
    ("cubic_c1 eval_s2_d1_1", -0.23978951315076263),
    ("cubic_c1 eval_s2_d1_2", 5.1567435060016047),
    ("cubic_c1 eval_s2_d2_0", -54.794734547422017),
    ("cubic_c1 eval_s2_d2_1", -214.94809126096317),
    ("cubic_c1 eval_s2_d2_2", 98.60485916975891),
    ("cubic_c1 trim_s_of_mid", 0.46467537933458769),
    ("quintic length", 6.889477399380989),
    ("quintic length_tol", 6.88954512519655),
    ("quintic length_ab", 3.4611750802239705),
    ("quintic length_ab_tol", 3.4611776726800887),
    ("quintic cpnts_ab", 3.4611750802239705),
    ("quintic cpnts_ab_tol", 3.4611776726800887),
    ("quintic abscissa_fwd", 0.36914185828726259),
    ("quintic abscissa_back", 0.89100560703566445),
    ("quintic abscissa_mid", 0.9029299442058184),
    ("quintic abscissa_tol_fwd", 0.36914185898546475),
    ("quintic abscissa_tol_back", 0.8500526718588135),
    ("quintic abscissa_guess_tol", 0.67690470182179141),
    ("quintic abscissa_guess", 0.47312936951590256),
    ("quintic cpnts_abscissa", 0.27898159484966301),
    ("quintic curvlin_length", 6.8895451252104163),
    ("quintic s_of_u0", 0.14014841451973961),
    ("quintic s_of_u1", 0.43238393190213054),
    ("quintic s_of_u2", 0.81232489209065928),
    ("quintic u_of_s0", 0.0),
    ("quintic u_of_s1", 0.027772598755051964),
    ("quintic u_of_s2", 0.36914656075342184),
    ("quintic u_of_s3", 0.54015064082073561),
    ("quintic u_of_s4", 0.55376974469701057),
    ("quintic u_of_s5", 0.96927482540261611),
    ("quintic u_of_s6", 0.99999999999999989),
    ("quintic u_of_s7", 0.16117204493151749),
    ("quintic nb_intervals_c2", 1.0),
    ("quintic interval_c2_1", 0.0),
    ("quintic interval_c2_2", 1.0),
    ("quintic trim_first", 0.14999999999999999),
    ("quintic trim_last", 0.84999999999999998),
    ("quintic eval_s0_d0_0", 1.2344079340161982),
    ("quintic eval_s0_d0_1", 0.39879838932673239),
    ("quintic eval_s0_d0_2", 0.18027471281156418),
    ("quintic eval_s0_d1_0", 6.5378078663768235),
    ("quintic eval_s0_d1_1", -1.0230723234515777),
    ("quintic eval_s0_d1_2", 1.9173480006611501),
    ("quintic eval_s0_d2_0", -3.2730350172640943),
    ("quintic eval_s0_d2_1", -8.4264167383288839),
    ("quintic eval_s0_d2_2", 6.6642259660408794),
    ("quintic eval_s1_d0_0", 3.2007548735881679),
    ("quintic eval_s1_d0_1", 0.66669995564052553),
    ("quintic eval_s1_d0_2", 0.55917953846342117),
    ("quintic eval_s1_d1_0", 6.7955502752996519),
    ("quintic eval_s1_d1_1", 1.1166167590475227),
    ("quintic eval_s1_d1_2", -0.19873475186700118),
    ("quintic eval_s1_d2_0", 3.0921221576350231),
    ("quintic eval_s1_d2_1", -19.537018581202418),
    ("quintic eval_s1_d2_2", -4.0390056720116148),
    ("quintic eval_s2_d0_0", 5.1012474277116375),
    ("quintic eval_s2_d0_1", 0.57568952822728747),
    ("quintic eval_s2_d0_2", 1.1430727644363166),
    ("quintic eval_s2_d1_0", 5.1262623582258273),
    ("quintic eval_s2_d1_1", -1.5579530078354058),
    ("quintic eval_s2_d1_2", 4.3312871865444569),
    ("quintic eval_s2_d2_0", -7.8732232567094531),
    ("quintic eval_s2_d2_1", -14.261811643343984),
    ("quintic eval_s2_d2_2", 4.1883566917946844),
    ("quintic trim_s_of_mid", 0.47022556442375285),
    ("bezier length", 6.3667584073752916),
    ("bezier length_tol", 6.3666680127172688),
    ("bezier length_ab", 3.787665237541745),
    ("bezier length_ab_tol", 3.7876652543220501),
    ("bezier cpnts_ab", 3.787665237541745),
    ("bezier cpnts_ab_tol", 3.7876652543220501),
    ("bezier abscissa_fwd", 0.33380766456954336),
    ("bezier abscissa_back", 0.81901319047664556),
    ("bezier abscissa_mid", 0.86220013858275923),
    ("bezier abscissa_tol_fwd", 0.33380766456951672),
    ("bezier abscissa_tol_back", 0.75954978588467925),
    ("bezier abscissa_guess_tol", 0.58586249732961049),
    ("bezier abscissa_guess", 0.4205697666414468),
    ("bezier cpnts_abscissa", 0.26028307654587607),
    ("bezier curvlin_length", 6.366668012619888),
    ("bezier s_of_u0", 0.12900406202340631),
    ("bezier s_of_u1", 0.47689059363228692),
    ("bezier s_of_u2", 0.87482110962582149),
    ("bezier u_of_s0", 0.0),
    ("bezier u_of_s1", 0.036765259741642346),
    ("bezier u_of_s2", 0.33380203702262351),
    ("bezier u_of_s3", 0.47534994038962619),
    ("bezier u_of_s4", 0.4863331581652181),
    ("bezier u_of_s5", 0.94921078062007991),
    ("bezier u_of_s6", 0.99999999999999989),
    ("bezier u_of_s7", 0.16266515926237338),
    ("bezier nb_intervals_c2", 1.0),
    ("bezier interval_c2_1", 0.0),
    ("bezier interval_c2_2", 1.0),
    ("bezier trim_first", 0.14999999999999999),
    ("bezier trim_last", 0.84999999999999998),
    ("bezier eval_s0_d0_0", 0.77708784767852934),
    ("bezier eval_s0_d0_1", 0.98730247148317407),
    ("bezier eval_s0_d0_2", 0.11271140591689242),
    ("bezier eval_s0_d1_0", 4.9538888240193684),
    ("bezier eval_s0_d1_1", 3.8590655317745806),
    ("bezier eval_s0_d1_2", 1.0493142158775153),
    ("bezier eval_s0_d2_0", 10.420014272327641),
    ("bezier eval_s0_d2_1", -14.378842074828125),
    ("bezier eval_s0_d2_2", 3.6874575121974305),
    ("bezier eval_s1_d0_0", 2.5510442976173473),
    ("bezier eval_s1_d0_1", 1.346592840951482),
    ("bezier eval_s1_d0_2", 0.47529298204586889),
    ("bezier eval_s1_d1_0", 6.1154139502747693),
    ("bezier eval_s1_d1_1", -1.3835899534775669),
    ("bezier eval_s1_d1_2", 1.1053744344510414),
    ("bezier eval_s1_d2_0", -2.7704248447671929),
    ("bezier eval_s1_d2_1", -12.691877713552438),
    ("bezier eval_s1_d2_2", -0.55914062429608258),
    ("bezier eval_s2_d0_0", 4.2555722182372335),
    ("bezier eval_s2_d0_1", 0.75245855718694932),
    ("bezier eval_s2_d0_2", 1.0317432514534848),
    ("bezier eval_s2_d1_0", 5.0443959537412972),
    ("bezier eval_s2_d1_1", -0.7804441602727844),
    ("bezier eval_s2_d1_2", 3.8051856666254551),
    ("bezier eval_s2_d2_0", -10.696741251686467),
    ("bezier eval_s2_d2_1", 23.121241989293296),
    ("bezier eval_s2_d2_2", 18.922450278458381),
    ("bezier trim_s_of_mid", 0.52243164405746689),
    ("circle length", 4.5),
    ("circle length_tol", 4.5),
    ("circle length_ab", 2.8800000000000003),
    ("circle length_ab_tol", 2.8800000000000003),
    ("circle cpnts_ab", 2.8800000000000003),
    ("circle cpnts_ab_tol", 2.8800000000000003),
    ("circle abscissa_fwd", 0.96599999999999997),
    ("circle abscissa_back", 1.7400000000000002),
    ("circle abscissa_mid", 1.7400000000000002),
    ("circle abscissa_tol_fwd", 0.96599999999999997),
    ("circle abscissa_tol_back", 1.6500000000000001),
    ("circle abscissa_guess_tol", 1.3799999999999999),
    ("circle abscissa_guess", 1.1099999999999999),
    ("circle cpnts_abscissa", 0.84000000000000008),
    ("circle curvlin_length", 4.5),
    ("circle s_of_u0", 0.099999999999999992),
    ("circle s_of_u1", 0.45000000000000007),
    ("circle s_of_u2", 0.90000000000000013),
    ("circle u_of_s0", 0.29999999999999999),
    ("circle u_of_s1", 0.39000000000000001),
    ("circle u_of_s2", 0.96599999999999975),
    ("circle u_of_s3", 1.1999999999999997),
    ("circle u_of_s4", 1.2179999999999997),
    ("circle u_of_s5", 1.9739999999999998),
    ("circle u_of_s6", 2.0999999999999996),
    ("circle u_of_s7", 0.65999999999999992),
    ("circle nb_intervals_c2", 1.0),
    ("circle interval_c2_1", 0.0),
    ("circle interval_c2_2", 1.0),
    ("circle trim_first", 0.14999999999999999),
    ("circle trim_last", 0.84999999999999998),
    ("circle eval_s0_d0_0", 2.9749805787434127),
    ("circle eval_s0_d0_1", 3.5327921299335845),
    ("circle eval_s0_d0_2", 3.0),
    ("circle eval_s0_d1_0", -2.7590258338804516),
    ("circle eval_s0_d1_1", 3.5549650417381433),
    ("circle eval_s0_d1_2", 0.0),
    ("circle eval_s0_d2_0", -6.3989370751286581),
    ("circle eval_s0_d2_1", -4.9662465009848127),
    ("circle eval_s0_d2_2", 0.0),
    ("circle eval_s1_d0_0", 1.9058943861916846),
    ("circle eval_s1_d0_1", 4.3300977149180655),
    ("circle eval_s1_d0_2", 3.0),
    ("circle eval_s1_d1_0", -4.1941758868525181),
    ("circle eval_s1_d1_1", 1.6306098951450323),
    ("circle eval_s1_d1_2", 0.0),
    ("circle eval_s1_d2_0", -2.9350978112610582),
    ("circle eval_s1_d2_1", -7.5495165963345325),
    ("circle eval_s1_d2_2", 0.0),
    ("circle eval_s2_d0_0", 0.57900638012730798),
    ("circle eval_s2_d0_1", 4.4642979470888839),
    ("circle eval_s2_d0_2", 3.0),
    ("circle eval_s2_d1_0", -4.4357363047599909),
    ("circle eval_s2_d1_1", -0.75778851577084561),
    ("circle eval_s2_d1_2", 0.0),
    ("circle eval_s2_d2_0", 1.3640193283875222),
    ("circle eval_s2_d2_1", -7.9843253485679835),
    ("circle eval_s2_d2_2", -0.0),
    ("circle trim_s_of_mid", 0.5),
    ("line length", 5.0),
    ("line length_tol", 5.0),
    ("line length_ab", 3.2000000000000002),
    ("line length_ab_tol", 3.2000000000000002),
    ("line cpnts_ab", 3.2000000000000002),
    ("line cpnts_ab_tol", 3.2000000000000002),
    ("line abscissa_fwd", 0.85000000000000009),
    ("line abscissa_back", 3.0),
    ("line abscissa_mid", 3.0),
    ("line abscissa_tol_fwd", 0.85000000000000009),
    ("line abscissa_tol_back", 2.75),
    ("line abscissa_guess_tol", 2.0),
    ("line abscissa_guess", 1.25),
    ("line cpnts_abscissa", 0.5),
    ("line curvlin_length", 5.0),
    ("line s_of_u0", 0.10000000000000001),
    ("line s_of_u1", 0.45000000000000001),
    ("line s_of_u2", 0.90000000000000002),
    ("line u_of_s0", -1.0),
    ("line u_of_s1", -0.75),
    ("line u_of_s2", 0.85000000000000009),
    ("line u_of_s3", 1.5),
    ("line u_of_s4", 1.55),
    ("line u_of_s5", 3.6500000000000004),
    ("line u_of_s6", 4.0),
    ("line u_of_s7", 0.0),
    ("line nb_intervals_c2", 1.0),
    ("line interval_c2_1", 0.0),
    ("line interval_c2_2", 1.0),
    ("line trim_first", 0.14999999999999999),
    ("line trim_last", 0.84999999999999998),
    ("line eval_s0_d0_0", 1.0),
    ("line eval_s0_d0_1", 1.0),
    ("line eval_s0_d0_2", 1.0),
    ("line eval_s0_d1_0", 1.6666666666666665),
    ("line eval_s0_d1_1", 3.333333333333333),
    ("line eval_s0_d1_2", 3.333333333333333),
    ("line eval_s0_d2_0", 0.0),
    ("line eval_s0_d2_1", 0.0),
    ("line eval_s0_d2_2", 0.0),
    ("line eval_s1_d0_0", 1.5),
    ("line eval_s1_d0_1", 2.0),
    ("line eval_s1_d0_2", 2.0),
    ("line eval_s1_d1_0", 1.6666666666666665),
    ("line eval_s1_d1_1", 3.333333333333333),
    ("line eval_s1_d1_2", 3.333333333333333),
    ("line eval_s1_d2_0", 0.0),
    ("line eval_s1_d2_1", 0.0),
    ("line eval_s1_d2_2", 0.0),
    ("line eval_s2_d0_0", 2.0),
    ("line eval_s2_d0_1", 3.0),
    ("line eval_s2_d0_2", 3.0),
    ("line eval_s2_d1_0", 1.6666666666666665),
    ("line eval_s2_d1_1", 3.333333333333333),
    ("line eval_s2_d1_2", 3.333333333333333),
    ("line eval_s2_d2_0", 0.0),
    ("line eval_s2_d2_1", 0.0),
    ("line eval_s2_d2_2", 0.0),
    ("line trim_s_of_mid", 0.5),
    ("ellipse length", 5.6771127037748421),
    ("ellipse length_tol", 5.6771127174640377),
    ("ellipse length_ab", 3.9194244943166909),
    ("ellipse length_ab_tol", 3.9194244943174001),
    ("ellipse cpnts_ab", 3.9194244943166909),
    ("ellipse cpnts_ab_tol", 3.9194244943174001),
    ("ellipse abscissa_fwd", 1.1957040548262752),
    ("ellipse abscissa_back", 2.0300845128325671),
    ("ellipse abscissa_mid", 1.8830419258041493),
    ("ellipse abscissa_tol_fwd", 1.1957040548269269),
    ("ellipse abscissa_tol_back", 1.9284603932330209),
    ("ellipse abscissa_guess_tol", 1.6383032764377201),
    ("ellipse abscissa_guess", 1.3529625608031037),
    ("ellipse cpnts_abscissa", 1.0513200877418947),
    ("ellipse curvlin_length", 5.6771127174519105),
    ("ellipse s_of_u0", 0.059853714313946718),
    ("ellipse s_of_u1", 0.38967986715589592),
    ("ellipse s_of_u2", 0.90948330286835022),
    ("ellipse u_of_s0", 0.20000000000000001),
    ("ellipse u_of_s1", 0.39568231426551698),
    ("ellipse u_of_s2", 1.1957040566177344),
    ("ellipse u_of_s3", 1.4487690243149205),
    ("ellipse u_of_s4", 1.4677939702023595),
    ("ellipse u_of_s5", 2.3188500960100789),
    ("ellipse u_of_s6", 2.5),
    ("ellipse u_of_s7", 0.82588776638034822),
    ("ellipse nb_intervals_c2", 1.0),
    ("ellipse interval_c2_1", 0.0),
    ("ellipse interval_c2_2", 1.0),
    ("ellipse trim_first", 0.14999999999999999),
    ("ellipse trim_last", 0.84999999999999998),
    ("ellipse eval_s0_d0_0", 2.033713774465248),
    ("ellipse eval_s0_d0_1", 0.88217987132333187),
    ("ellipse eval_s0_d0_2", 0.0),
    ("ellipse eval_s0_d1_0", -5.3263338904321742),
    ("ellipse eval_s0_d1_1", 1.9646312871090774),
    ("ellipse eval_s0_d1_2", 0.0),
    ("ellipse eval_s0_d2_0", -3.0911652751910683),
    ("ellipse eval_s0_d2_1", -8.380492805041623),
    ("ellipse eval_s0_d2_2", 0.0),
    ("ellipse eval_s1_d0_0", 0.36517404994199543),
    ("ellipse eval_s1_d0_1", 1.1910766835598092),
    ("ellipse eval_s1_d0_2", 0.0),
    ("ellipse eval_s1_d1_0", -5.6702944360894039),
    ("ellipse eval_s1_d1_1", 0.27815430017845422),
    ("ellipse eval_s1_d1_2", 0.0),
    ("ellipse eval_s1_d2_0", -0.21454128719667098),
    ("ellipse eval_s1_d2_1", -4.3735159453665524),
    ("ellipse eval_s1_d2_2", 0.0),
    ("ellipse eval_s2_d0_0", -1.3299305129297403),
    ("ellipse eval_s2_d0_1", 1.0756419352761168),
    ("ellipse eval_s2_d0_2", 0.0),
    ("ellipse eval_s2_d1_0", -5.5691842427445426),
    ("ellipse eval_s2_d1_1", -1.1017239568150081),
    ("ellipse eval_s2_d1_2", 0.0),
    ("ellipse eval_s2_d2_0", 1.0931300032606996),
    ("ellipse eval_s2_d2_1", -5.5257420443413272),
    ("ellipse eval_s2_d2_2", 0.0),
    ("ellipse trim_s_of_mid", 0.44846590881008042),
    ("periodic length", 8.884299408517716),
    ("periodic length_tol", 8.8842994098358581),
    ("periodic length_ab", 5.6620365468603886),
    ("periodic length_ab_tol", 5.6620365468657088),
    ("periodic cpnts_ab", 5.6691626094079757),
    ("periodic cpnts_ab_tol", 5.6620365468657319),
    ("periodic abscissa_fwd", 1.8114563119048532),
    ("periodic abscissa_back", 3.980950162109866),
    ("periodic abscissa_mid", 3.990727018668498),
    ("periodic abscissa_tol_fwd", 1.8114563119874174),
    ("periodic abscissa_tol_back", 3.6841667202036694),
    ("periodic abscissa_guess_tol", 2.7657728083891464),
    ("periodic abscissa_guess", 2.1540259198812421),
    ("periodic cpnts_abscissa", 1.4918284254389866),
    ("periodic curvlin_length", 8.8842994098358581),
    ("periodic s_of_u0", 0.10861236614294381),
    ("periodic s_of_u1", 0.4736804475126497),
    ("periodic s_of_u2", 0.89887016876869297),
    ("periodic u_of_s0", 0.0),
    ("periodic u_of_s1", 0.23138995258601636),
    ("periodic u_of_s2", 1.8114563122309177),
    ("periodic u_of_s3", 2.3548724367182778),
    ("periodic u_of_s4", 2.3945088159309091),
    ("periodic u_of_s5", 4.6566671579335086),
    ("periodic u_of_s6", 4.9999999999999982),
    ("periodic u_of_s7", 0.96903490054140906),
    ("periodic nb_intervals_c2", 1.0),
    ("periodic interval_c2_1", 0.0),
    ("periodic interval_c2_2", 1.0),
    ("periodic trim_first", 0.14999999999999999),
    ("periodic trim_last", 0.84999999999999998),
    ("periodic eval_s0_d0_0", 2.5140541904842331),
    ("periodic eval_s0_d0_1", 1.786421110700694),
    ("periodic eval_s0_d0_2", 0.66572267452967382),
    ("periodic eval_s0_d1_0", -2.2884759169687077),
    ("periodic eval_s0_d1_1", 8.5777949040477406),
    ("periodic eval_s0_d1_2", 0.33924705622397999),
    ("periodic eval_s0_d2_0", -93.959853656751591),
    ("periodic eval_s0_d2_1", -22.704522237284479),
    ("periodic eval_s0_d2_2", -59.750340460375121),
    ("periodic eval_s1_d0_0", 0.31260045848358597),
    ("periodic eval_s1_d0_1", 2.163455888020831),
    ("periodic eval_s1_d0_2", 0.2371116156437553),
    ("periodic eval_s1_d1_0", -7.1821975122751489),
    ("periodic eval_s1_d1_1", -5.2059119372411251),
    ("periodic eval_s1_d1_2", 0.49527345978639437),
    ("periodic eval_s1_d2_0", 18.00940134949774),
    ("periodic eval_s1_d2_1", -23.937623718957695),
    ("periodic eval_s1_d2_2", 9.5501109715412618),
    ("periodic eval_s2_d0_0", 0.13827566211988784),
    ("periodic eval_s2_d0_1", 0.17637303385959913),
    ("periodic eval_s2_d0_2", 0.088184788650720353),
    ("periodic eval_s2_d1_0", 8.271858560156776),
    ("periodic eval_s2_d1_1", -2.8998674686627433),
    ("periodic eval_s2_d1_2", -1.4484131415319226),
    ("periodic eval_s2_d2_0", 16.032706924823731),
    ("periodic eval_s2_d2_1", 36.951920511519091),
    ("periodic eval_s2_d2_2", 17.581041687932668),
    ("periodic trim_s_of_mid", 0.53652460254470624),
    ("segment_bspline length", 5.0),
    ("segment_bspline length_tol", 5.0),
    ("segment_bspline length_ab", 3.2000000000000002),
    ("segment_bspline length_ab_tol", 3.2000000000000002),
    ("segment_bspline cpnts_ab", 3.2000000000000002),
    ("segment_bspline cpnts_ab_tol", 3.2000000000000002),
    ("segment_bspline abscissa_fwd", 0.73999999999999999),
    ("segment_bspline abscissa_back", 1.6000000000000001),
    ("segment_bspline abscissa_mid", 1.6000000000000001),
    ("segment_bspline abscissa_tol_fwd", 0.73999999999999999),
    ("segment_bspline abscissa_tol_back", 1.5),
    ("segment_bspline abscissa_guess_tol", 1.2),
    ("segment_bspline abscissa_guess", 0.90000000000000002),
    ("segment_bspline cpnts_abscissa", 0.59999999999999998),
    ("segment_bspline curvlin_length", 5.0),
    ("segment_bspline s_of_u0", 0.10000000000000001),
    ("segment_bspline s_of_u1", 0.45000000000000001),
    ("segment_bspline s_of_u2", 0.90000000000000002),
    ("segment_bspline u_of_s0", 0.0),
    ("segment_bspline u_of_s1", 0.10000000000000001),
    ("segment_bspline u_of_s2", 0.73999999999999999),
    ("segment_bspline u_of_s3", 1.0),
    ("segment_bspline u_of_s4", 1.02),
    ("segment_bspline u_of_s5", 1.8599999999999999),
    ("segment_bspline u_of_s6", 1.9999999999999998),
    ("segment_bspline u_of_s7", 0.40000000000000002),
    ("segment_bspline nb_intervals_c2", 1.0),
    ("segment_bspline interval_c2_1", 0.0),
    ("segment_bspline interval_c2_2", 1.0),
    ("segment_bspline trim_first", 0.14999999999999999),
    ("segment_bspline trim_last", 0.84999999999999998),
    ("segment_bspline eval_s0_d0_0", 1.6000000000000001),
    ("segment_bspline eval_s0_d0_1", 0.80000000000000004),
    ("segment_bspline eval_s0_d0_2", 2.0),
    ("segment_bspline eval_s0_d1_0", 3.0),
    ("segment_bspline eval_s0_d1_1", 4.0),
    ("segment_bspline eval_s0_d1_2", 0.0),
    ("segment_bspline eval_s0_d2_0", 0.0),
    ("segment_bspline eval_s0_d2_1", 0.0),
    ("segment_bspline eval_s0_d2_2", 0.0),
    ("segment_bspline eval_s1_d0_0", 2.5),
    ("segment_bspline eval_s1_d0_1", 2.0),
    ("segment_bspline eval_s1_d0_2", 2.0),
    ("segment_bspline eval_s1_d1_0", 3.0),
    ("segment_bspline eval_s1_d1_1", 4.0),
    ("segment_bspline eval_s1_d1_2", 0.0),
    ("segment_bspline eval_s1_d2_0", 0.0),
    ("segment_bspline eval_s1_d2_1", 0.0),
    ("segment_bspline eval_s1_d2_2", 0.0),
    ("segment_bspline eval_s2_d0_0", 3.3999999999999999),
    ("segment_bspline eval_s2_d0_1", 3.1999999999999997),
    ("segment_bspline eval_s2_d0_2", 2.0),
    ("segment_bspline eval_s2_d1_0", 3.0),
    ("segment_bspline eval_s2_d1_1", 4.0),
    ("segment_bspline eval_s2_d1_2", 0.0),
    ("segment_bspline eval_s2_d2_0", 0.0),
    ("segment_bspline eval_s2_d2_1", 0.0),
    ("segment_bspline eval_s2_d2_2", 0.0),
    ("segment_bspline trim_s_of_mid", 0.5),
    ("segment_bezier length", 3.4871191548325395),
    ("segment_bezier length_tol", 3.4871191548325395),
    ("segment_bezier length_ab", 2.2317562590928253),
    ("segment_bezier length_ab_tol", 2.2317562590928253),
    ("segment_bezier cpnts_ab", 2.2317562590928253),
    ("segment_bezier cpnts_ab_tol", 2.2317562590928253),
    ("segment_bezier abscissa_fwd", 0.39600000000000002),
    ("segment_bezier abscissa_back", 0.73999999999999999),
    ("segment_bezier abscissa_mid", 0.73999999999999999),
    ("segment_bezier abscissa_tol_fwd", 0.39600000000000002),
    ("segment_bezier abscissa_tol_back", 0.69999999999999996),
    ("segment_bezier abscissa_guess_tol", 0.58000000000000007),
    ("segment_bezier abscissa_guess", 0.46000000000000008),
    ("segment_bezier cpnts_abscissa", 0.34000000000000002),
    ("segment_bezier curvlin_length", 3.4871191548325395),
    ("segment_bezier s_of_u0", 0.10000000000000001),
    ("segment_bezier s_of_u1", 0.45000000000000007),
    ("segment_bezier s_of_u2", 0.90000000000000002),
    ("segment_bezier u_of_s0", 0.10000000000000001),
    ("segment_bezier u_of_s1", 0.14000000000000001),
    ("segment_bezier u_of_s2", 0.39600000000000002),
    ("segment_bezier u_of_s3", 0.5),
    ("segment_bezier u_of_s4", 0.5079999999999999),
    ("segment_bezier u_of_s5", 0.84399999999999986),
    ("segment_bezier u_of_s6", 0.8999999999999998),
    ("segment_bezier u_of_s7", 0.26000000000000001),
    ("segment_bezier nb_intervals_c2", 1.0),
    ("segment_bezier interval_c2_1", 0.0),
    ("segment_bezier interval_c2_2", 1.0),
    ("segment_bezier trim_first", 0.14999999999999999),
    ("segment_bezier trim_last", 0.84999999999999998),
    ("segment_bezier eval_s0_d0_0", -0.21999999999999997),
    ("segment_bezier eval_s0_d0_1", 0.26000000000000001),
    ("segment_bezier eval_s0_d0_2", 0.78000000000000003),
    ("segment_bezier eval_s0_d1_0", 2.4000000000000004),
    ("segment_bezier eval_s0_d1_1", 0.80000000000000004),
    ("segment_bezier eval_s0_d1_2", 2.4000000000000004),
    ("segment_bezier eval_s0_d2_0", 0.0),
    ("segment_bezier eval_s0_d2_1", 0.0),
    ("segment_bezier eval_s0_d2_2", 0.0),
    ("segment_bezier eval_s1_d0_0", 0.5),
    ("segment_bezier eval_s1_d0_1", 0.5),
    ("segment_bezier eval_s1_d0_2", 1.5),
    ("segment_bezier eval_s1_d1_0", 2.4000000000000004),
    ("segment_bezier eval_s1_d1_1", 0.80000000000000004),
    ("segment_bezier eval_s1_d1_2", 2.4000000000000004),
    ("segment_bezier eval_s1_d2_0", 0.0),
    ("segment_bezier eval_s1_d2_1", 0.0),
    ("segment_bezier eval_s1_d2_2", 0.0),
    ("segment_bezier eval_s2_d0_0", 1.2199999999999998),
    ("segment_bezier eval_s2_d0_1", 0.73999999999999988),
    ("segment_bezier eval_s2_d0_2", 2.2199999999999998),
    ("segment_bezier eval_s2_d1_0", 2.4000000000000004),
    ("segment_bezier eval_s2_d1_1", 0.80000000000000004),
    ("segment_bezier eval_s2_d1_2", 2.4000000000000004),
    ("segment_bezier eval_s2_d2_0", 0.0),
    ("segment_bezier eval_s2_d2_1", 0.0),
    ("segment_bezier eval_s2_d2_2", 0.0),
    ("segment_bezier trim_s_of_mid", 0.5),
];

fn p(x: f64, y: f64, z: f64) -> Point {
    Point::new(x, y, z)
}

/// The curves of the OCCT driver, with the ranges their adaptors take.
fn curves() -> Vec<(&'static str, Curve, f64, f64)> {
    let frame = |o: Point| Frame::new(o, Direction::Z, Direction::X);
    vec![
        (
            "cubic_c1",
            Curve::BSpline(BSplineCurve::new(
                vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(2.0, -1.0, 1.0), p(3.0, 3.0, 0.0), p(4.0, 0.0, 2.0), p(5.0, 1.0, -1.0), p(6.0, 2.0, 0.0), p(7.0, 0.0, 1.0)],
                vec![0.0, 1.0, 2.0, 3.0, 4.0],
                vec![4, 1, 2, 1, 4],
                3,
            )),
            0.0,
            4.0,
        ),
        (
            "quintic",
            Curve::BSpline(BSplineCurve::new(
                vec![p(0.0, 0.0, 0.0), p(1.0, 1.0, 0.0), p(2.0, -1.0, 0.5), p(3.0, 2.0, 1.0), p(4.0, 0.0, 0.0), p(5.0, 1.0, 1.0), p(6.0, 0.0, 2.0)],
                vec![0.0, 0.5, 1.0],
                vec![6, 1, 6],
                5,
            )),
            0.0,
            1.0,
        ),
        ("bezier", Curve::Bezier(BezierCurve::new(vec![p(0.0, 0.0, 0.0), p(1.0, 2.0, 0.0), p(3.0, 2.0, 1.0), p(4.0, 0.0, 0.0), p(5.0, 1.0, 2.0)])), 0.0, 1.0),
        ("circle", Curve::Circle(Circle::new(frame(p(1.0, 2.0, 3.0)), 2.5)), 0.3, 2.1),
        ("line", Curve::Line(Line::new(p(1.0, 1.0, 1.0), Direction::new(1.0, 2.0, 2.0))), -1.0, 4.0),
        ("ellipse", Curve::Ellipse(Ellipse::new(frame(Point::ORIGIN), 3.0, 1.2)), 0.2, 2.5),
        (
            "periodic",
            Curve::BSpline(BSplineCurve::new_periodic(
                vec![p(0.0, 0.0, 0.0), p(2.0, 0.0, 0.0), p(3.0, 2.0, 1.0), p(1.0, 3.0, 0.0), p(-1.0, 1.0, 0.5)],
                vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                vec![1, 1, 1, 1, 1, 1],
                3,
            )),
            0.0,
            5.0,
        ),
        ("segment_bspline", Curve::BSpline(BSplineCurve::new(vec![p(1.0, 0.0, 2.0), p(4.0, 4.0, 2.0)], vec![0.0, 2.0], vec![2, 2], 1)), 0.0, 2.0),
        ("segment_bezier", Curve::Bezier(BezierCurve::new(vec![p(-1.0, 0.0, 0.0), p(2.0, 1.0, 3.0)])), 0.1, 0.9),
    ]
}

/// The driver's calls, in its order.
fn run(name: &str, curve: &Curve, f: f64, l: f64, out: &mut Vec<(String, f64)>) {
    let mut put = |key: &str, v: f64| out.push((format!("{name} {key}"), v));
    let c = GeomAdaptor::new(curve, f, l);
    let r = l - f;
    let (a, b) = (f + 0.13 * r, f + 0.77 * r);
    let len = length(&c, None);
    put("length", len);
    put("length_tol", length(&c, Some(1e-6)));
    put("length_ab", length_between(&c, a, b, None));
    put("length_ab_tol", length_between(&c, a, b, Some(1e-8)));
    put("cpnts_ab", integrated_length(&c, a, b, None));
    put("cpnts_ab_tol", integrated_length(&c, a, b, Some(1e-8)));
    put("abscissa_fwd", AbscissaPoint::new(&c, 0.37 * len, f).unwrap().parameter().unwrap());
    put("abscissa_back", AbscissaPoint::new(&c, -0.2 * len, l).unwrap().parameter().unwrap());
    put("abscissa_mid", AbscissaPoint::new(&c, 0.5 * len, f + 0.3 * r).unwrap().parameter().unwrap());
    put("abscissa_tol_fwd", AbscissaPoint::with_tolerance(1e-7, &c, 0.37 * len, f).parameter().unwrap());
    put("abscissa_tol_back", AbscissaPoint::with_tolerance(1e-7, &c, -0.25 * len, l).parameter().unwrap());
    put("abscissa_guess_tol", AbscissaPoint::with_guess_and_tolerance(&c, 0.6 * len, f, f + 0.5 * r, 1e-9).parameter().unwrap());
    put("abscissa_guess", AbscissaPoint::with_guess(&c, 0.45 * len, f, f + 0.2 * r).parameter().unwrap());
    put("cpnts_abscissa", CurveAbscissaPoint::new(&c, 0.3 * len, f, 1e-9).parameter().unwrap());

    let mut function = CurvlinFunc::new(&c, 1e-8);
    put("curvlin_length", function.get_length());
    for (i, u) in [0.1, 0.45, 0.9].into_iter().enumerate() {
        put(&format!("s_of_u{i}"), function.get_s_parameter(f + u * r));
    }
    for (i, s) in [0.0, 0.05, 0.37, 0.5, 0.51, 0.93, 1.0, 0.2].into_iter().enumerate() {
        put(&format!("u_of_s{i}"), function.get_u_parameter_on(&c, s).unwrap());
    }
    let n = function.nb_intervals(Continuity::C2);
    put("nb_intervals_c2", n as f64);
    for (i, s) in function.intervals(Continuity::C2).into_iter().enumerate() {
        put(&format!("interval_c2_{}", i + 1), s);
    }
    function.trim(0.15, 0.85, 1e-9).unwrap();
    put("trim_first", function.first_parameter());
    put("trim_last", function.last_parameter());
    for (i, s) in [0.2, 0.5, 0.8].into_iter().enumerate() {
        for order in 0..=2 {
            let v = function.eval_case1(s, order).unwrap().unwrap();
            for (j, x) in v.into_iter().enumerate() {
                put(&format!("eval_s{i}_d{order}_{j}"), x);
            }
        }
    }
    put("trim_s_of_mid", function.get_s_parameter(f + 0.5 * r));
}

/// Every value against OCCT's, to 1e-12 relative (absolute below 1).
#[test]
fn agrees_with_occt() {
    let mut values = Vec::new();
    for (name, curve, f, l) in curves() {
        run(name, &curve, f, l, &mut values);
    }
    assert_eq!(values.len(), OCCT_VALUES.len());
    let mut worst = (0.0, String::new());
    for ((key, value), (occt_key, occt)) in values.iter().zip(OCCT_VALUES) {
        assert_eq!(key, occt_key);
        let deviation = (value - occt).abs() / occt.abs().max(1.0);
        if deviation > worst.0 {
            worst = (deviation, key.clone());
        }
        assert!(deviation <= 1e-12, "{key}: {value:.17e} against OCCT's {occt:.17e} ({deviation:.1e})");
    }
    eprintln!("largest deviation from OCCT: {:.1e} ({})", worst.0, worst.1);
}

/// GetUParameter undoes GetSParameter, and both run from end to end.
#[test]
fn arc_length_round_trip() {
    for (name, curve, f, l) in curves() {
        let c = GeomAdaptor::new(&curve, f, l);
        let function = CurvlinFunc::new(&c, 1e-10);
        assert!((function.get_s_parameter(l) - 1.0).abs() < 1e-12, "{name}");
        for k in 0..=20 {
            let s = k as f64 / 20.0;
            let u = function.get_u_parameter(s).unwrap();
            assert!((function.get_s_parameter(u) - s).abs() < 1e-8, "{name} at s = {s}");
        }
    }
}

/// Arc lengths found by the inverse come back as lengths.
#[test]
fn abscissa_inverts_length() {
    for (name, curve, f, l) in curves() {
        let c = GeomAdaptor::new(&curve, f, l);
        let len = length(&c, Some(1e-10));
        for k in 1..10 {
            let target = len * k as f64 / 10.0;
            let u = AbscissaPoint::with_tolerance(1e-10, &c, target, f).parameter().unwrap();
            assert!((length_between(&c, f, u, Some(1e-10)) - target).abs() < 1e-7 * len, "{name}: {k}/10");
        }
    }
}

/// The inputs OCCT refuses.
#[test]
fn errors_where_occt_raises() {
    let curve = Curve::Line(Line::new(Point::ORIGIN, Direction::X));
    let c = GeomAdaptor::new(&curve, 0.0, 1e-9);
    assert!(AbscissaPoint::new(&c, 0.5, 0.0).is_err());
    let c = GeomAdaptor::new(&curve, 0.0, 1.0);
    let mut function = CurvlinFunc::new(&c, 1e-8);
    assert!(function.get_u_parameter(1.5).is_err());
    assert!(function.trim(-0.1, 0.5, 1e-9).is_err());
}

/// Gauss quadrature with n points integrates polynomials of degree 2n - 1 exactly.
#[test]
fn gauss_quadrature_is_exact_for_polynomials() {
    for order in 1..=24usize {
        let degree = 2 * order as i32 - 1;
        let f = |x: f64| x.powi(degree) + 1.0;
        let exact = (2f64.powi(degree + 1) - (-1f64).powi(degree + 1)) / (degree + 1) as f64 + 3.0;
        let value = gauss_single_integration(&f, -1.0, 2.0, order);
        assert!((value - exact).abs() < 1e-12 * exact.abs().max(1.0) * 2f64.powi(degree), "order {order}: {value} against {exact}");
    }
}
