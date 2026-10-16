# DRAW conversion report

OCCT: `V8_0_1-42-g3d097a0328`; groups: blend, boolean, bugs, chamfer, cr, draft, evolved, feat, fillet2d, geometry, heal, helix, lowalgos, mesh, mkface, nproject, offset, pipe, prism, prj, sewing, thrusection

- scripts considered: 16420
- converted: 10245 (2858 ignored for missing data), 6373 of them as YAML cases (cases/occt/draw)
- skipped: 6175
- distinct DRAW commands used by converted scripts: 342

## Skip reasons

| reason | scripts | examples |
|---|---|---|
| tcl:if | 1683 | boolean/bsection/N4, boolean/cells_test/D5, boolean/cells_test/J3 |
| grid end file: tcl:regexp | 1395 | heal/direct_faces/A1, heal/direct_faces/A2, heal/direct_faces/A3 |
| no commands | 536 | bugs/moddata_1/bug22554_1, bugs/moddata_1/bug22554_2, bugs/moddata_1/bug22554_3 |
| tcl:source | 328 | boolean/gdml_private/A1, boolean/gdml_private/A2, boolean/gdml_private/A3 |
| grid begin file | 289 | blend/tolblend_simple/A1, blend/tolblend_simple/A2, blend/tolblend_simple/A3 |
| tcl:foreach | 288 | blend/complex/A6, boolean/bcut_complex/F2, boolean/bsection/M8 |
| tcl:list | 284 | bugs/modalg_6/bug26647, bugs/modalg_6/bug28795, bugs/modalg_7/bug32863 |
| tcl:regexp | 269 | boolean/bsection/M3, bugs/fclasses/bug23774, bugs/fclasses/bug24137 |
| tcl:lrange | 244 | bugs/modalg_6/bug21507, bugs/moddata_3/bug28182, bugs/moddata_3/bug28183 |
| tcl:catch | 180 | bugs/fclasses/bug23192_1, bugs/fclasses/bug23192_2, bugs/fclasses/bug28829 |
| tcl:llength | 61 | bugs/modalg_1/buc60663_2, bugs/modalg_1/buc60708, bugs/modalg_1/bug13140 |
| tcl:for | 60 | bugs/fclasses/bug24863_1, bugs/fclasses/bug7287_1, bugs/fclasses/bug7287_2 |
| data-exchange:stepread | 46 | boolean/removefeatures/A1, boolean/removefeatures/A2, boolean/removefeatures/A3 |
| tcl:lindex | 45 | boolean/bfuse_simple/L1, bugs/fclasses/bug1488, bugs/modalg_2/bug22833 |
| proc:CheckIsFeatureRemoved | 43 | boolean/removefeatures/A9, boolean/removefeatures/B1, boolean/removefeatures/B2 |
| tcl:lsearch | 43 | thrusection/specific/A1, thrusection/specific/A2, thrusection/specific/A3 |
| tcl:string | 42 | bugs/modalg_2/bug22428, bugs/modalg_2/bug23089, bugs/modalg_5/bug25124_6 |
| dynamic set | 42 | bugs/modalg_2/bug415, bugs/modalg_7/bug30396_1, lowalgos/2dinter/bug27976 |
| tcl:eval | 41 | bugs/modalg_6/bug24097, bugs/modalg_6/bug26701, bugs/modalg_6/bug27283 |
| non-constant expr | 18 | bugs/modalg_2/bug442, bugs/modalg_6/bug27537, bugs/moddata_1/buc60852 |
| grid end file: tcl:if | 17 | heal/drop_small_edges/A1, heal/drop_small_edges/A2, heal/drop_small_edges/A3 |
| proc:perform_offset_increasing | 17 | offset/shape_type_i_c/XS4, offset/shape_type_i_c/XS5, offset/shape_type_i_c/XS6 |
| tcl:proc | 14 | bugs/modalg_2/bug427_1, bugs/modalg_2/bug427_2, bugs/modalg_2/bug427_3 |
| unknown variable | 13 | boolean/simplify/A3, boolean/simplify/A4, bugs/fclasses/bug1281 |
| viewer:vdump | 11 | bugs/fclasses/bug30804, bugs/modalg_7/bug28144, bugs/moddata_3/bug26359 |
| grid end file: unknown variable | 11 | bugs/splitshape/case1, bugs/splitshape/case10, bugs/splitshape/case11 |
| grid end file: tcl:proc | 11 | heal/drop_small_solids/M1, heal/drop_small_solids/M2, heal/drop_small_solids/M3 |
| proc:CheckPerform | 11 | heal/checkshape/bug27814_1, heal/checkshape/bug27814_10, heal/checkshape/bug27814_11 |
| {*} expansion | 10 | boolean/history/A7, bugs/modalg_6/bug26180, bugs/modalg_6/bug27878_1 |
| tcl:incr | 9 | boolean/bsection/N5, boolean/bsection/N6, boolean/bsection/N7 |
| viewer:vhlr | 9 | bugs/modalg_5/bug24842_axo, bugs/modalg_5/bug24842_back, bugs/modalg_5/bug24842_bottom |
| tcl:split | 9 | bugs/moddata_1/bug22489_1, bugs/moddata_1/bug22489_2, bugs/moddata_2/bug23248 |
| proc:val3d | 9 | geometry/beziecurve/B4, geometry/beziecurve/B5, geometry/bsplinecurve/C2 |
| data-exchange:testreadstep | 7 | bugs/modalg_5/bug25477_2, bugs/modalg_6/bug28028, bugs/modalg_7/bug30196 |
| proc:val2d | 7 | geometry/2dbeziecurve/B1, geometry/2dbsplinecurve/B7, geometry/circle/C2 |
| data-exchange:ReadStep | 6 | bugs/modalg_5/bug24140, bugs/modalg_5/bug24746, bugs/modalg_6/bug26308 |
| data-exchange:igesread | 5 | bugs/modalg_4/bug8370, bugs/modalg_5/bug24012, bugs/modalg_5/bug24074 |
| array variable | 5 | bugs/modalg_7/bug30273, bugs/moddata_3/bug31136_1, bugs/moddata_3/bug31136_2 |
| data-exchange:param | 5 | bugs/mesh/bug25281, bugs/mesh/bug25287, bugs/mesh/bug25837_1 |
| grid end file: array variable | 5 | heal/elementary_to_revolution/A1, heal/elementary_to_revolution/A2, heal/elementary_to_revolution/A3 |
| proc:perform_offset_multi_with_ref | 5 | offset/shape_type_i_c_multi/C2, offset/shape_type_i_c_multi/C3, offset/shape_type_i_c_multi/XI1 |
| tcl:while | 4 | bugs/modalg_6/bug26193, bugs/modalg_6/bug26208, bugs/heal/bug25333 |
| viewer:vreadpixel | 4 | bugs/modalg_6/bug26305_1, bugs/modalg_6/bug26305_2, bugs/modalg_6/bug26305_3 |
| tcl:open | 4 | bugs/modalg_7/bug29939, bugs/modalg_7/bug30269, bugs/modalg_7/optimal_bndbox_exception |
| checks only in the group's end script | 4 | offset/faces_type_i/C9, offset/shape_type_i/C5, offset/bugs/bug26030 |
| tcl:encoding | 3 | bugs/fclasses/bug25367_brep, bugs/fclasses/bug25367_igs, bugs/fclasses/bug63 |
| viewer:vmoveto | 2 | bugs/modalg_2/bug21909, bugs/mesh/bug25469_3 |
| viewer:vselect | 2 | bugs/modalg_2/bug22781_2, bugs/modalg_2/bug22781_4 |
| data-exchange:testreadiges | 2 | bugs/modalg_5/bug25410, bugs/heal/bug27078 |
| proc:checkarea | 2 | bugs/modalg_6/bug26607, bugs/heal/bug26644 |
| data-exchange:igesbrep | 2 | bugs/moddata_1/bug15570, bugs/heal/bug26671 |
| proc:CheckExtResult | 2 | lowalgos/extcc/bug29465_1, lowalgos/extcc/bug29465_2 |
| draw_deviations.json: Known to fail: the volume maker's union of these faces fails after about 45 s (a section edge left without a pcurve), and the test then runs about 90 s in all; not converted until the faces' union succeeds. | 1 | boolean/volumemaker/G1 |
| draw_deviations.json: checks OCCT's shape sharing, not geometry (nbshapes counts a vertex and its moved copy once) | 1 | bugs/modalg_2/bug22588 |
| data-exchange:stepwrite | 1 | bugs/modalg_6/bug22454 |
| dynamic data file name | 1 | bugs/modalg_6/bug26717 |
| viewer:vhlrtype | 1 | bugs/modalg_7/bug28780 |
| viewer:vcomputehlr | 1 | bugs/modalg_7/bug29405 |
| data-exchange:readstl | 1 | bugs/modalg_7/bug29524 |
| tcl:unset | 1 | bugs/moddata_3/bug25481 |
| data-exchange:newmodel | 1 | bugs/moddata_3/bug25971 |
| data-exchange:brepiges | 1 | bugs/heal/bug23722 |
| tcl:concat | 1 | heal/reshape/A1 |

## Commands by use

| command | uses | scripts |
|---|---|---|
| `restore` | 9669 | 6957 |
| `checkprops` | 8074 | 7435 |
| `dset` | 5768 | 4077 |
| `checkshape` | 3425 | 2814 |
| `checknbshapes` | 3392 | 3010 |
| `bop` | 2746 | 2685 |
| `box` | 2522 | 1670 |
| `explode` | 1751 | 1439 |
| `trotate` | 1486 | 1110 |
| `ttranslate` | 1310 | 910 |
| `pcylinder` | 1157 | 948 |
| `checksection` | 1089 | 1079 |
| `vertex` | 1088 | 262 |
| `edge` | 872 | 169 |
| `offsetparameter` | 852 | 852 |
| `offsetload` | 852 | 852 |
| `offsetperform` | 852 | 852 |
| `checkmaxtol` | 842 | 836 |
| `mkface` | 824 | 347 |
| `sewing` | 752 | 747 |
| `bcut` | 739 | 587 |
| `checkfreebounds` | 739 | 738 |
| `wire` | 726 | 498 |
| `bopcommon` | 709 | 706 |
| `checkfaults` | 700 | 700 |
| `DT_SplitAngle` | 697 | 697 |
| `plane` | 695 | 263 |
| `bopfuse` | 635 | 600 |
| `bfuse` | 632 | 540 |
| `bopcut` | 630 | 617 |
| `boptuc` | 622 | 622 |
| `mkoffset_fraction` | 621 | 621 |
| `offset_counts` | 621 | 621 |
| `mkplane` | 571 | 422 |
| `copy` | 522 | 167 |
| `prism` | 519 | 349 |
| `mkedge` | 483 | 302 |
| `generated` | 445 | 42 |
| `pcone` | 420 | 329 |
| `profile` | 406 | 290 |
| `blend` | 354 | 337 |
| `nexplode` | 328 | 306 |
| `tscale` | 320 | 290 |
| `polyline` | 302 | 249 |
| `bcommon` | 299 | 294 |
| `psphere` | 298 | 239 |
| `thrusections` | 292 | 290 |
| `offsetonface` | 290 | 49 |
| `unifysamedom` | 270 | 256 |
| `bclearobjects` | 269 | 255 |
| `baddobjects` | 269 | 255 |
| `bfillds` | 269 | 255 |
| `bcleartools` | 263 | 249 |
| `circle` | 238 | 169 |
| `featperform` | 220 | 218 |
| `bsection` | 219 | 211 |
| `checktrinfo` | 186 | 173 |
| `tclean` | 182 | 168 |
| `trim` | 179 | 105 |
| `addsweep` | 169 | 108 |
| `incmesh` | 169 | 153 |
| `checklength` | 165 | 158 |
| `bopsection` | 152 | 151 |
| `pipe` | 149 | 137 |
| `fixshape` | 146 | 133 |
| `settolerance` | 145 | 45 |
| `featprism` | 143 | 141 |
| `baddtools` | 142 | 133 |
| `depouille` | 131 | 131 |
| `offset` | 128 | 87 |
| `bcadd` | 120 | 72 |
| `nurbsconvert` | 119 | 113 |
| `renamevar` | 116 | 107 |
| `ellipse` | 112 | 108 |
| `bbuild` | 110 | 104 |
| `ptorus` | 108 | 93 |
| `mksweep` | 108 | 108 |
| `buildsweep` | 108 | 108 |
| `mkevol` | 103 | 103 |
| `updatevol` | 103 | 103 |
| `buildevol` | 103 | 103 |
| `beziercurve` | 101 | 93 |
| `revol` | 100 | 94 |
| `mkvolume` | 98 | 98 |
| `add` | 97 | 65 |
| `bbop` | 95 | 72 |
| `featperformval` | 93 | 92 |
| `setflags` | 89 | 51 |
| `mksurface` | 89 | 75 |
| `compound` | 87 | 67 |
| `line` | 85 | 43 |
| `prj` | 81 | 81 |
| `tcopy` | 80 | 42 |
| `featrevol` | 80 | 80 |
| `bcbuild` | 78 | 78 |
| `cylinder` | 78 | 52 |
| `offsetshapesimple` | 78 | 77 |
| `bcremoveall` | 74 | 59 |
| `bsplinecurve` | 72 | 71 |
| `dval` | 65 | 22 |
| `checkreal` | 59 | 20 |
| `DT_ClosedSplit` | 59 | 59 |
| `DT_ToBspl` | 58 | 58 |
| `fsameparameter` | 55 | 49 |
| `thickshell` | 53 | 53 |
| `parabola` | 52 | 51 |
| `extsurf` | 52 | 52 |
| `triangles` | 50 | 50 |
| `orientation` | 49 | 31 |
| `point` | 48 | 11 |
| `project` | 47 | 46 |
| `savehistory` | 47 | 43 |
| `cone` | 46 | 35 |
| `checkdump` | 46 | 46 |
| `featlf` | 41 | 40 |
| `hyperbola` | 41 | 40 |
| `shape` | 40 | 32 |
| `convert` | 37 | 34 |
| `featdprism` | 36 | 36 |
| `sphere` | 35 | 27 |
| `mkcurve` | 35 | 24 |
| `bopcheck` | 34 | 25 |
| `updatetolerance` | 33 | 17 |
| `openoffset` | 32 | 31 |
| `setsweep` | 31 | 31 |
| `revsurf` | 30 | 30 |
| `halfspace` | 26 | 23 |
| `bfuzzyvalue` | 26 | 21 |
| `DT_SplitByNumber` | 26 | 26 |
| `trimv` | 25 | 12 |
| `binrestore` | 25 | 17 |
| `to3d` | 24 | 4 |
| `sameparameter` | 23 | 21 |
| `tolerance` | 23 | 16 |
| `polyvertex` | 22 | 19 |
| `wedge` | 21 | 21 |
| `distmini` | 21 | 16 |
| `nproject` | 21 | 21 |
| `bopargcheck` | 20 | 16 |
| `translate` | 20 | 19 |
| `bounding` | 19 | 15 |
| `mkoffset` | 19 | 17 |
| `pcurve` | 19 | 17 |
| `beziersurf` | 19 | 19 |
| `bfuseblend` | 18 | 18 |
| `nbshapes` | 18 | 14 |
| `brestore` | 18 | 14 |
| `chamf` | 17 | 8 |
| `btolx` | 16 | 10 |
| `reverse` | 16 | 16 |
| `middlepath` | 15 | 15 |
| `maxtolerance` | 14 | 9 |
| `bcremove` | 14 | 10 |
| `2dtranslate` | 14 | 14 |
| `bcaddall` | 13 | 13 |
| `bsplit` | 13 | 13 |
| `splitshape` | 13 | 13 |
| `2dbsplinecurve` | 13 | 13 |
| `v2d` | 13 | 13 |
| `pbsplinecurve` | 13 | 13 |
| `2dbeziercurve` | 13 | 10 |
| `featrf` | 13 | 13 |
| `GetShape` | 12 | 6 |
| `smirror` | 12 | 12 |
| `fixwgaps` | 12 | 12 |
| `encoderegularity` | 11 | 11 |
| `bcmakecontainers` | 10 | 10 |
| `bcremoveint` | 10 | 10 |
| `bapibop` | 10 | 5 |
| `mat` | 10 | 7 |
| `extrema` | 9 | 8 |
| `offsetshape` | 9 | 9 |
| `modified` | 9 | 5 |
| `addslide` | 9 | 9 |
| `fillet2d` | 9 | 9 |
| `lmirror` | 9 | 8 |
| `invert` | 8 | 8 |
| `result` | 8 | 5 |
| `evolved` | 8 | 8 |
| `rotate` | 8 | 7 |
| `trinfo` | 8 | 6 |
| `tricheck` | 8 | 5 |
| `pscale` | 8 | 8 |
| `topoload` | 7 | 7 |
| `fillet` | 7 | 7 |
| `concatC0wire` | 7 | 4 |
| `svalue` | 7 | 3 |
| `mk2dcurve` | 7 | 5 |
| `offsetcompshape` | 7 | 7 |
| `bsplinesurf` | 7 | 7 |
| `l` | 7 | 7 |
| `2dlmirror` | 7 | 7 |
| `side` | 6 | 3 |
| `drawcont` | 6 | 6 |
| `baddcompound` | 6 | 5 |
| `Open` | 6 | 6 |
| `bounds` | 6 | 5 |
| `removeinternals` | 6 | 4 |
| `2dpmirror` | 6 | 6 |
| `2drotate` | 6 | 6 |
| `pmirror` | 6 | 6 |
| `tmirror` | 5 | 5 |
| `dall` | 5 | 5 |
| `emptycopy` | 5 | 5 |
| `uiso` | 5 | 5 |
| `bsimplify` | 5 | 3 |
| `scalexyz` | 5 | 4 |
| `reflectlines` | 5 | 5 |
| `houtl` | 5 | 4 |
| `hfill` | 5 | 4 |
| `hload` | 5 | 4 |
| `2dcvalue` | 5 | 3 |
| `xdistcs` | 5 | 4 |
| `OCC30391` | 5 | 5 |
| `BUC60652` | 5 | 5 |
| `torus` | 5 | 5 |
| `proj` | 5 | 5 |
| `2dpscale` | 5 | 5 |
| `qcircle` | 5 | 3 |
| `breducetolerance` | 4 | 3 |
| `bopaddpcs` | 4 | 2 |
| `computetolerance` | 4 | 3 |
| `brollingball` | 4 | 4 |
| `hprj` | 4 | 4 |
| `hsetprj` | 4 | 4 |
| `hupdate` | 4 | 4 |
| `hhide` | 4 | 4 |
| `hres2d` | 4 | 4 |
| `baddctools` | 4 | 4 |
| `bnondestructive` | 4 | 4 |
| `boporigin` | 4 | 1 |
| `DT_ShapeDivide` | 4 | 4 |
| `tuyau` | 4 | 3 |
| `ssolid` | 4 | 1 |
| `incudeg` | 4 | 4 |
| `incvdeg` | 4 | 4 |
| `movep` | 4 | 4 |
| `viso` | 3 | 3 |
| `fixsmallfaces` | 3 | 2 |
| `appro` | 3 | 1 |
| `appsurf` | 3 | 1 |
| `circ2d3Tan` | 3 | 3 |
| `treverse` | 2 | 2 |
| `hlr` | 2 | 1 |
| `build3d` | 2 | 2 |
| `coordload` | 2 | 2 |
| `subshape` | 2 | 2 |
| `scaleshape` | 2 | 2 |
| `OCC606` | 2 | 2 |
| `offset2dcurve` | 2 | 2 |
| `offsetcurve` | 2 | 2 |
| `bclassify` | 2 | 2 |
| `parameters` | 2 | 1 |
| `firsthole` | 2 | 2 |
| `intersect` | 2 | 2 |
| `brunparallel` | 2 | 2 |
| `bglue` | 2 | 2 |
| `numshapes` | 2 | 2 |
| `BUC60843` | 2 | 2 |
| `segment` | 2 | 2 |
| `tocanon` | 2 | 2 |
| `bclear` | 2 | 2 |
| `OCC26930` | 2 | 2 |
| `fillcurves` | 2 | 1 |
| `insertpole` | 2 | 2 |
| `insertknot` | 2 | 2 |
| `qline` | 2 | 2 |
| `polygonprops` | 2 | 1 |
| `2dintersect` | 2 | 2 |
| `cirtang` | 2 | 2 |
| `bcutblend` | 1 | 1 |
| `bapisplit` | 1 | 1 |
| `BUC60944` | 1 | 1 |
| `OCC24834` | 1 | 1 |
| `OCC26485` | 1 | 1 |
| `OCC23429` | 1 | 1 |
| `OCC369` | 1 | 1 |
| `OCC527` | 1 | 1 |
| `OCC363` | 1 | 1 |
| `solidorientation` | 1 | 1 |
| `OCC23952sweep` | 1 | 1 |
| `OCC23952intersect` | 1 | 1 |
| `OCC24667` | 1 | 1 |
| `connectedges` | 1 | 1 |
| `concatwire` | 1 | 1 |
| `cprj` | 1 | 1 |
| `proximity` | 1 | 1 |
| `checkcurveonsurf` | 1 | 1 |
| `lprops` | 1 | 1 |
| `cvalue` | 1 | 1 |
| `bhaspc` | 1 | 1 |
| `OCC26270` | 1 | 1 |
| `OCC26313` | 1 | 1 |
| `bapibuild` | 1 | 1 |
| `OCC26396` | 1 | 1 |
| `OCC26446` | 1 | 1 |
| `OCC26553` | 1 | 1 |
| `hlrin3d` | 1 | 1 |
| `correctnormals` | 1 | 1 |
| `mkpoint` | 1 | 1 |
| `samerange` | 1 | 1 |
| `fuseedge` | 1 | 1 |
| `hremove` | 1 | 1 |
| `u` | 1 | 1 |
| `OCC28784` | 1 | 1 |
| `OCC32744` | 1 | 1 |
| `bopbsolid` | 1 | 1 |
| `OCC566` | 1 | 1 |
| `BUC60792` | 1 | 1 |
| `OCC14376` | 1 | 1 |
| `projonplane` | 1 | 1 |
| `OCC367` | 1 | 1 |
| `fixsmalledges` | 1 | 1 |
| `BUC60811` | 1 | 1 |
| `draft` | 1 | 1 |
| `OCC23945` | 1 | 1 |
| `wexplo` | 1 | 1 |
| `trimu` | 1 | 1 |
| `setunotperiodic` | 1 | 1 |
| `2dextrema` | 1 | 1 |
| `splitc12d` | 1 | 1 |
| `approxcurveonsurf` | 1 | 1 |
| `segsur` | 1 | 1 |
| `exchuv` | 1 | 1 |
| `bsmooth` | 1 | 1 |
| `tessellate` | 1 | 1 |
| `FAILED` | 1 | 1 |
| `removeloc` | 1 | 1 |
| `fixsmall` | 1 | 1 |
| `OCC7570` | 1 | 1 |
| `gplate` | 1 | 1 |
| `chamfer2d` | 1 | 1 |
| `incdeg` | 1 | 1 |
| `remknot` | 1 | 1 |
| `polygon2d` | 1 | 1 |
| `addpolygonnode` | 1 | 1 |
| `2dapprox` | 1 | 1 |
| `2dinterpole` | 1 | 1 |
| `disp` | 1 | 1 |
| `mkoffsetshape` | 1 | 1 |
| `clpoles` | 1 | 1 |
| `clknots` | 1 | 1 |
