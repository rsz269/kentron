# GTest porting status

Hand translations of OCCT GTests live in `tests/gtest_<toolkit>/`, one module per C++ file.
In scope: computational geometry and topology toolkits. Out of scope: TKernel (OS/collections/strings),
visualization, data exchange, OCAF, TKHLR, TKExpress.

Translation rules: tests use kentron's own API, not OCCT names; `docs/occt-mapping.md` maps
one to the other (extend it as new concepts are ported). Each module keeps the name of its
OCCT source file and each OCCT test keeps a counterpart, even when the translation is
trivial, so coverage stays traceable.

- `Handle(Geom_X)` becomes a plain owned value; `DownCast` becomes a `match` on the
  `Curve` / `Surface` enums; a null handle becomes `Option`.
- Topology operates on a `topology::Model` created in the test and passed to every call.
- `EXPECT_THROW` on a constructor becomes `assert!(...try_new(...).is_none())`; other
  exceptions become panics checked with `assert_panics!`.
- Tests relying on mutating an object through one handle and observing the change through
  another have no direct equivalent: restructure them and say so in a comment.

## Deviations from OCCT

Tests whose expectation intentionally differs from OCCT's:

| test | OCCT expects | kentron expects | why |
|---|---|---|---|
| `gtest_tkbrep::brep_tools_gtrsf_modification::tolerance_tracks_current_transformation` (shear case), `surface_located_bezier_transforms_evaluated_geometry` | tolerance unchanged (0.1); 0.2 x 4 (largest diagonal entry) | 0.1 x golden ratio; 0.2 x largest singular value | a tolerance must grow by the largest stretch (max singular value) to cover the deformed tolerance zone; OCCT uses the spectral radius |
| `gtest_tkprim::brep_prim_api_make_torus::partial_torus_{zero_to180, minus180_to_zero, 180_to360, positive_range, negative_range}`, `root_cause_height_ordering_assumption` | invalid shape | valid shape | OCCT asserts a documented limitation of its implementation; the solid is well defined |

| `gtest_tkbrep::brep_graph_compact::history_mode_records_mapping`, `history_survives_compact_with_remapped_ids` | `Compact:Remap` history events | no events; ids in the history are rewritten, `CompactReport::renumbered` counts renumbered nodes | compaction is bookkeeping, not a modeling operation |
| `gtest_tkbrep::brep_graph_deduplicate::canonicalize_surfaces_rewrites_without_history`, `history_records_merge_events` and other dedup history counts | one `Dedup:Canonicalize*` event per shared surface/curve, `Dedup:Merge*` per merged node | no event for geometry sharing; one `"deduplicate"` `Replaced` event per merged node | sharing geometry does not change any topology node |
| `gtest_tkbrep::brep_graph_deduplicate::idempotent_*` | second run may report rewrites again | second run reports none | already-shared geometry is not rewritten |
| `gtest_tkbrep::brep_graph_sparse_model::compact_remove_part_shape_empties_part`, `remove_*` | removal leaves dangling references until `CleanupRemovedReferences`; removing a part's solid retires its occurrence | `remove_node` cleans up immediately; the part stays, with no shape | no half-removed state; a part's shape is not an instance |
| `gtest_tkg3d::geom_api_extrema_curve_curve` | `Distance(i)` in solver order | extrema sorted by increasing distance | a stable, meaningful order |
| `gtest_tkg2d::geom2d_api_inter_curve_curve::intersecting_ellipses_points_lie_on_both` | roots checked with OCCT's Newton solver on a trigonometric equation | each point checked on both ellipses | same property, without porting `math_NewtonFunctionRoot` |
| `gtest_tkgeombase::extrema_pc2d_conics::circle::multiple_periods_return_every_extremum` (and every `PointCurveSearch` result) | extrema in solver order (period by period, minima first) | extrema sorted by increasing parameter | a stable, meaningful order |
| `gtest_tkbo` (every `validate_result`) | surface area within ±5000 | within 1e-4 relative | a ±5000 tolerance accepts nearly any result, and hid the two expectation bugs below |
| `gtest_tkbo::brep_algo_api_cut::complex_profile_*` (H4, H5, H8, H9, I1–I4, I7–I9) | profiles built without the `p 0 0 -1 1 0 0` their comments quote | `p` restored | the expected areas are those of the DRAW scripts, which have it (H4: 98000 with it, 95000 without) |
| `gtest_tkbo::brep_algo_api_fuse::nurbs_box_plus_narrow_box_c4`, `c5`, `b7`, `b8`, `c8`, `c9` | areas 6.75, 7.75, 6.25 | 7, 10.5, 8.5 | the tests build other boxes than the DRAW scripts they cite but kept those scripts' areas; these are the areas of the boxes built |
| `gtest_tkbo::brep_algo_api_fuse::large_prism_with_oblong_prism_right_side_j5`, `j6`, `j7` | area 180.283 | 178 | 180.283 is the DRAW script's, fusing an oblong prism with rounded ends; the test fuses a 1 x 2 x 2 box instead |
| `gtest_tkbrep::brep_graph_revisions::adding_an_instance_bumps_the_assembly` | child product's revision bumped too | only the assembly's | a product does not store its instances |
| `gtest_tkmath::bvh_box::area_degenerate_box2_d` | area of a flat rectangle: width + height (5) | 0 | OCCT's `Area` serves its surface area heuristic; kentron's `BoundingBox2::area` is the area |
| `gtest_tkg2d::geom2d_offset_curve::rational_bspline_offset_approximation_keeps_length` | length 3.1415850499856011 within 1e-6 | π within 1e-5 | the expected value carries OCCT's particular approximation error (7.6e-6 short of π, within its 1e-4 tolerance); another approximation within tolerance has another length. The check against π still catches the bug (a wrong offset) |
| `gtest_tkgeombase::project::cone_circle_parallel_axes_gives_line`, `cone_circle_opposite_normal_gives_reversed_line` | a coaxial circle of radius 2 (off the cone) maps to a `v = const` line | the circles are put on the cone | `analytic_pcurve` is for curves lying on the surface and rejects others (callers such as `unify` rely on that); OCCT's `ProjLib_Cone` projects any coaxial circle |

## Not ported

OCCT tests that check implementation details with no counterpart in kentron. Each test
module's header comment lists the ones it skips.

| test(s) | reason |
|---|---|
| `TopoDS_Builder_Test`: `TShapeFlags`, `Add_SetsModifiedFlag`, `Remove_SetsModifiedFlag`; `TopoDS_TShape_Test`: `FlagSettersGetters`, `Orientable_DifferentTypes`, `FlagsIndependent`; `TopoDS_Edge_Test`: `BUC60828_InfiniteFlag` | OCCT `TShape` bookkeeping flags (`Free`, `Modified`, `Checked`, `Infinite`, ...) |
| `TopoDS_Iterator_Test`: `OCC30708_1_InitializeWithNullShape`, `ValueThrowsWhenNotMore` | kentron has no null shapes; Rust iterators return `None` past the end |
| `BRepTools_GTrsfModificationTest`: `Curve2dAndVertexParameterPreserveParametrization`, `PolygonOnTriangulation_PreservesDeflection`; `BRepTools_TrsfModificationTest`: `ParameterlessPolygonOnTriangulationDoesNotThrow` | pcurves / mesh polylines bound to a face without a wire, vertex parameters on edges: partial topology OCCT allows and kentron's coedge model does not |
| binary round trips in `BRepTools_Test` | OCCT's binary BRep format (BinTools) |
| BRepGraph_EventBus, Lock, PermissionUpdate, MutGuard, DeferredInvalidation, CacheRegistry, LayerIterator, LayerTopoSupplement, ItemId, TypedIdDispatch | BRepGraph framework mechanics (event dispatch, locking, RAII mutation guards, plug-in layers, generic id unions) |
| BRepGraph_Reconstruct, BRepGraphInc_PopulateRegression; import-pipeline tests in BRepGraph_Build and elsewhere | conversion between BRepGraph and OCCT's TopoDS shapes, which kentron does not have |
| `BRepGraph_LayerHistory_Test`: `Absorb_*` (6) | bridge from OCCT's `BRepTools_History` keyed by `TopoDS_Shape` |
| `BRepGraph_Compact_Test`: `FindNodeStillWorksAfterCompact`, `CacheMesh_*` (3), `Compact_Preserves{DeletedItemUidHistory, ItemUidHistoryMappings}` | TopoDS bindings; the mesh cache layer; item-UID history (removal does not write history in kentron, and a `History<StableId>` is unaffected by compaction by construction) |
| `BRepGraph_Deduplicate_Test`: duplicate-pcurve tests (7), `RoundTrip_*`, `Build_SharedTFace_*`, `ParallelBuild_*`, `RestoresHistory*` (3) | coedges outside any wire; TopoDS round trips; parallel import; kentron never toggles the history's enabled flag |
| `BRepGraph_Copy_Test`: runtime mesh-cache and cache-service tests (4), `CopyBox_SharedGeometry`, `CopyGraph_PreservesSupplementAttachmentsAndUids`, `CopyGraph_NonEmptyTargetUsesSetterPathForPersistentMesh` | caches; geometry sharing between graphs (geometry is by value); OCCT's side layer for internal vertices; copy-into-non-empty is `Model::import` |
| `BRepGraph_Transform_Test`: `LocationOnly_NoCopyGeom`, `TransformNode_CopyGeomAndMesh_DropsRuntimeCache`, `PerformIntoNonEmptyTargetDoesNotMoveExistingVertices` | shared-geometry placement mode (`Shape::moved`); caches; no transform-into-target operation |
| `BRepGraph_Validate_Test`: `DocumentRootReferencedByOccurrence_Detected`, `Audit_DetectsSharedOwnedUse`, `RemoveNodeMaintainsActiveCount`, `Audit_ValidatesCoEdgeUIDs*`, `BoundsCheck_InvalidIndex` (merged into `orphan_curve_slot_is_an_error`) | root products are derived, not stored; shared slots are legal (deduplication); internal counters |
| `BRepGraph_Assembly_Test`: `Build_*AutoCreatesRootProduct`, `Add_*` (6), `ShapesView_*` (5), `LinkProducts_ParentOccurrenceMustMatchParentProduct`, `MutProduct_RAII`, `MutInvalidAssemblyDefs_*`, `NodeId_Helpers`, `UID_IsAssembly`, repeated count/relations variants | products are created explicitly, not on import; TopoDS reconstruction; occurrence-path disambiguation (kentron instances have no parent occurrence); mutation guards; id-kind helpers |
| BRepGraph_CacheMesh | the runtime mesh cache and its staleness tracking; kentron stores meshes on faces and coedges and exposes `Model::revision` for freshness checks |
| `BRepGraph_Polygon_Test`: `*_Roundtrip_*`, `Polygon3D_Captured_WhenPresent`, `VertexPointRepresentations_*`; `BRepGraph_EdgeCases_Test`: null-shape and rebuild tests, `ParallelBuild_*` | TopoDS round trips; vertex parameters on curves/surfaces (not stored in kentron); a rebuilt model is a new `Model` |
| `BRepGraph_Test`: `Reconstruct*`, `Shape_*`, `HasOriginalShape_*`, allocator, parallel-build and rebuild-UID tests, `Remove*_PrunesOrphan*`, `RemoveShapeRoot_*`, `AddPCurve_*`, `DetectToleranceConflicts_*` | TopoDS reconstruction; detaching never removes implicitly; coedges outside wires; vacuous on a box |
| `BRepGraphInc_Test`: `Populate_*`, `Reconstruct_*`, `*_RoundTrip_*` via TopoDS, `Storage_*`, `Relations_*` (reverse-index maintenance), supplement-layer routing, `ParityOrientation*` | TopoDS import/export; storage internals and the incremental reverse index; OCCT's side layer for internal/external sub-shapes |
| `BRepGraph_ShapesViewImport_Test`: `Mutable*Definition`, `MutInvalidTopologyDefs_*`, `AddWireToFace_*`, `AddWithHistory_*`, `SkipsFirstNode` and duplicates; `BRepGraph_ScenarioMatrix_Test`: reconstruct / populate round trips, `*_HasNoOriginal`, `*_ShapeReconstructsSubEdge` | mutation guards; the same edges in two wires of one face; TopoDS import and original-shape bindings |
| `BRepGraph_RefId_Test`: `RefUID_*`, `RefDomain_Stamp*`, `MutFaceRef_*`, `StaleRefUID_*` | no default/generation-tagged stable ids; stamps are for nodes only; a rebuilt model is a new `Model` |
| Geometry `Copy*` tests (copy constructors, `Copy()`, copy independence), `WeightsArray_*`, null-handle and 1-based-index checks | geometry values are plain `Clone` data; weights are `Option<&[f64]>` / `weight(i)` |
| Geom2d/Geom `*Eval` curve and surface families (AHT/T-Bezier, spirals, involute, sine wave, helix, quadric patches) | recent non-core OCCT additions; candidates for later |
| `*GridEval*`, `*Hash*`, `GeomAdaptor_Transformed*` | batch-evaluation, hashing and adaptor-caching internals |
| `BRepGraph_VersionStamp_Test`: `DefaultStamp_IsInvalid`, `IsStale_InvalidStamp_ReturnsTrue`, `IsStale_DeferredMode_TracksCorrectly`, `GenericItemUID_*`, `StampOf_RemovedUse_*`, `GraphGUID_*`, `ToGUID_*` | no default/invalid stamps (`Option<Stamp>`); deferred invalidation; stamps of use and geometry records; GUID persistence |
| `BRepGraph_MutationGen_Test`: `OwnGen_DeferredMode`, `SubtreeGen_Deferred*`, `EditorCreationRegistersOwnGeneration` | deferred invalidation; in kentron new records start at revision 0 (OCCT counts creation-time wiring as edits) |
| `GeomLib_Test`: `CancelDenominatorDerivative_*`, `ParameterArrays_*`, `DensifyArray_*`, `FuseIntervals_*`, `IsIsoLine_*`, `GTransformLine_*` | approximation-internal reparametrization and array helpers; `isIsoLine` is internal; kentron has no 2D affine transformation |
| `Hermit_Test`, `LProp_CurAndInf_Test`, AdvApp2Var_*, AppCont_ContMatrices | internals of rational approximation and of the AdvApp2Var/AppCont approximation engines; `LProp_CurAndInf` is a result container (kentron returns a sorted `Vec<CurvatureFeature>`) |
| `IntAna_IntQuadQuad_Test`: `NextCurveBoundaryConditions` | out-of-range index exceptions; curve links are `Option<(usize, bool)>` |
| GeomBndLib_*: `*_CompareWithBndLib`, `AdaptorConstructor_*`, `TrimmedCurveHandle_*`, `*_UsesBasisSpecialization`, `SamplingHelpers*` | compare OCCT's two bounding-box implementations (or its adaptor/handle entry points) with each other; kentron has one, `geometry::bounds` |
| ExtremaPC/ExtremaPC2d: `TransformedCurve_CompoundScale_*PreservesParameter`, `ScaledTransformedOffsetUsesParameterSafeFallback`; dispatcher-vs-specialist comparisons, null-handle, ownership and result-storage tests; `*RefreshesPartitionAfterMutation`, `AnalyticalDispatchCopiesPrimitiveGeometry`, `GridUtilities*`; ExtremaPC2d_DistanceFunction, ExtremaPC2d_GridEvaluator, ExtremaPC2d.hxx helpers | OCCT's transformed-curve adaptor (keeps the untransformed parameter); one search instead of per-type evaluators; evaluator caches and internal helpers |
| ExtremaPC2d_OtherCurve | needs `Geom2dEval_SineWaveCurve` (Eval families are not ported) |
| BVH_BuildQueue, BVH_QuickSorter, BVH_RadixSorter, BVH_DistanceField | internals of OCCT's tree builders (parallel build queue, sorting by Morton codes) and its voxel distance field |
| BVH_Ray | a ray type caching the inverse direction; kentron's rays are an origin and a direction (`BoundingBox::ray_interval`, `Bvh::along_ray`, tested with BVH_Tools) |
| BVH_Triangulation | OCCT's primitive set (vertices, element indices, reordering during builds); a `Bvh` is built over boxes, a triangle's being `BoundingBox::from_points` |
| C++ id mechanics in BRepGraph_NodeId/RefId/RepId (`++`, hashing across types, invalid enum values) | no Rust counterpart |

## TKMath (120 files)

| file | tests | lines | status |
|---|---|---|---|
| Bnd_B2_Test.cxx | 18 | 264 |  |
| Bnd_B3_Test.cxx | 21 | 338 |  |
| Bnd_BoundSortBox_Test.cxx | 9 | 331 |  |
| Bnd_Box2d_Test.cxx | 38 | 443 |  |
| Bnd_Box_Test.cxx | 47 | 881 |  |
| Bnd_OBB_Test.cxx | 6 | 106 |  |
| Bnd_Range_Test.cxx | 54 | 511 |  |
| Bnd_Sphere_Test.cxx | 20 | 224 |  |
| BSplCLib_Cache_Test.cxx | 7 | 522 |  |
| BSplCLib_Test.cxx | 5 | 76 |  |
| BSplSLib_Cache_Test.cxx | 8 | 978 |  |
| BSplSLib_Test.cxx | 5 | 93 |  |
| BVH_BinnedBuilder_Test.cxx | 14 | 298 | done (2 not ported) |
| BVH_Box_Test.cxx | 86 | 1134 | done (15 not ported: 4D, single precision, constexpr, corner setters) |
| BVH_BuildQueue_Test.cxx | 17 | 389 | not ported |
| BVH_DistanceField_Test.cxx | 1 | 133 | not ported |
| BVH_LinearBuilder_Test.cxx | 15 | 330 | done (3 not ported) |
| BVH_PairDistance_Test.cxx | 2 | 240 | done |
| BVH_QuickSorter_Test.cxx | 13 | 322 | not ported |
| BVH_RadixSorter_Test.cxx | 35 | 449 | not ported |
| BVH_Ray_Test.cxx | 19 | 268 | not ported |
| BVH_ShapeSelect_Test.cxx | 2 | 244 | done |
| BVH_SpatialMedianBuilder_Test.cxx | 15 | 398 | done (3 not ported) |
| BVH_SweepPlaneBuilder_Test.cxx | 10 | 314 | done |
| BVH_Tools_Test.cxx | 45 | 776 | done (5 single-precision not ported) |
| BVH_Traverse_Test.cxx | 21 | 666 | done (2 null-tree tests not ported) |
| BVH_Tree_Test.cxx | 23 | 399 | done (22 not ported: node-by-node construction) |
| BVH_Triangulation_Test.cxx | 16 | 328 | not ported |
| Convert_CircleToBSplineCurve_Test.cxx | 6 | 163 |  |
| Convert_CompBezierCurvesToBSplineCurve_Test.cxx | 7 | 197 |  |
| Convert_CompPolynomialToPoles_Test.cxx | 4 | 153 |  |
| Convert_ConeToBSplineSurface_Test.cxx | 4 | 112 |  |
| Convert_CylinderToBSplineSurface_Test.cxx | 5 | 132 |  |
| Convert_EllipseToBSplineCurve_Test.cxx | 4 | 139 |  |
| Convert_GridPolynomialToPoles_Test.cxx | 2 | 121 |  |
| Convert_HyperbolaToBSplineCurve_Test.cxx | 3 | 139 |  |
| Convert_ParabolaToBSplineCurve_Test.cxx | 3 | 100 |  |
| Convert_SphereToBSplineSurface_Test.cxx | 6 | 141 |  |
| Convert_TorusToBSplineSurface_Test.cxx | 7 | 148 |  |
| CSLib_Test.cxx | 36 | 883 |  |
| ElCLib_Test.cxx | 12 | 414 |  |
| gp_Ax3_Test.cxx | 1 | 130 | done |
| gp_Circ_Test.cxx | 8 | 96 | done |
| gp_Dir2d_Test.cxx | 14 | 138 | done |
| gp_Dir_Test.cxx | 18 | 198 | done |
| gp_GTrsf_Test.cxx | 2 | 55 | done |
| gp_Lin_Test.cxx | 13 | 137 | done |
| gp_Mat_Test.cxx | 1 | 35 | done |
| gp_Pln_Test.cxx | 11 | 124 | done |
| gp_Pnt2d_Test.cxx | 14 | 171 | done |
| gp_Pnt_Test.cxx | 16 | 181 | done |
| gp_Quaternion_Test.cxx | 5 | 186 | done |
| gp_Torus_Test.cxx | 1 | 118 | done |
| gp_Trsf_Test.cxx | 1 | 39 | done |
| gp_Vec2d_Test.cxx | 17 | 180 | done |
| gp_Vec_Test.cxx | 18 | 187 | done |
| gp_XYZ_Test.cxx | 17 | 160 | done |
| math_BFGS_Test.cxx | 13 | 561 |  |
| math_BissecNewton_Test.cxx | 12 | 298 |  |
| math_BracketedRoot_Test.cxx | 14 | 244 |  |
| math_BracketMinimum_Test.cxx | 13 | 317 |  |
| math_BrentMinimum_Test.cxx | 16 | 320 |  |
| math_ComputeKronrodPointsAndWeights_Test.cxx | 1 | 36 |  |
| math_Crout_Test.cxx | 11 | 340 |  |
| math_DirectPolynomialRoots_Test.cxx | 19 | 466 |  |
| math_DoubleTab_Test.cxx | 11 | 215 |  |
| math_EigenValuesSearcher_Test.cxx | 30 | 1228 |  |
| math_FRPR_Test.cxx | 14 | 483 |  |
| math_FunctionAllRoots_Test.cxx | 16 | 417 |  |
| math_FunctionRoots_Test.cxx | 15 | 482 |  |
| math_FunctionRoot_Test.cxx | 12 | 396 |  |
| math_FunctionSetRoot_Test.cxx | 17 | 637 |  |
| math_GaussLeastSquare_Test.cxx | 14 | 444 |  |
| math_Gauss_Test.cxx | 12 | 414 |  |
| math_GlobOptMin_Test.cxx | 16 | 656 |  |
| math_Householder_Test.cxx | 10 | 330 |  |
| MathInteg_Comparison_Test.cxx | 21 | 628 |  |
| MathInteg_New_Test.cxx | 24 | 507 |  |
| math_Integration_Test.cxx | 16 | 440 |  |
| MathInteg_Test.cxx | 60 | 1171 |  |
| math_Jacobi_Test.cxx | 11 | 367 |  |
| MathLin_Comparison_Test.cxx | 12 | 517 |  |
| MathLin_EigenSearch_Test.cxx | 7 | 353 |  |
| MathLin_Test.cxx | 50 | 1624 |  |
| math_Matrix_Test.cxx | 11 | 649 |  |
| math_NewtonFunctionRoot_Test.cxx | 16 | 414 |  |
| math_NewtonFunctionSetRoot_Test.cxx | 14 | 513 |  |
| math_NewtonMinimum_Test.cxx | 15 | 528 |  |
| MathOpt_1D_Comparison_Test.cxx | 13 | 444 |  |
| MathOpt_1D_Test.cxx | 27 | 451 |  |
| MathOpt_FRPR_Test.cxx | 12 | 569 |  |
| MathOpt_GlobOpt_Test.cxx | 24 | 756 |  |
| MathOpt_NDim_Comparison_Test.cxx | 12 | 816 |  |
| MathOpt_NDim_Test.cxx | 25 | 931 |  |
| MathOpt_Newton_Test.cxx | 15 | 776 |  |
| MathOpt_PSO_Test.cxx | 28 | 928 |  |
| MathPoly_Comparison_Test.cxx | 18 | 412 |  |
| MathPoly_Laguerre_Test.cxx | 20 | 426 |  |
| MathPoly_Test.cxx | 47 | 608 |  |
| math_Powell_Test.cxx | 11 | 405 |  |
| math_PSO_Test.cxx | 13 | 531 |  |
| MathRoot_Comparison_Test.cxx | 15 | 611 |  |
| MathRoot_Multiple_Test.cxx | 27 | 606 |  |
| MathRoot_Test.cxx | 51 | 987 |  |
| MathRoot_Trig_Test.cxx | 50 | 715 |  |
| math_SVD_Test.cxx | 13 | 461 |  |
| MathSys_Comparison_Test.cxx | 24 | 1184 |  |
| MathSys_LM_Test.cxx | 30 | 1091 |  |
| MathSys_Newton2D_Test.cxx | 12 | 474 |  |
| MathSys_Newton3D_Test.cxx | 10 | 393 |  |
| MathSys_Newton4D_Test.cxx | 8 | 340 |  |
| math_TrigonometricFunctionRoots_Test.cxx | 16 | 298 |  |
| MathUtils_Bracket_Test.cxx | 16 | 364 |  |
| MathUtils_Functor_Test.cxx | 49 | 756 |  |
| math_Uzawa_Test.cxx | 14 | 411 |  |
| math_Vector_Test.cxx | 31 | 905 |  |
| PLib_HermitJacobi_Test.cxx | 9 | 294 |  |
| PLib_JacobiPolynomial_Test.cxx | 11 | 366 |  |
| PLib_Test.cxx | 19 | 528 |  |
| Poly_ArrayOfNodes_Test.cxx | 2 | 61 |  |

## TKG2d (33 files)

| file | tests | lines | status |
|---|---|---|---|
| Adaptor2d_Line2d_Test.cxx | 23 | 227 | merged into geom2d_line (adaptor layer not ported) |
| Adaptor2d_OffsetCurve_Test.cxx | 20 | 211 | merged into geom2d_offset_curve |
| Geom2dAdaptor_Curve_Test.cxx | 14 | 239 | merged into geom2d_trimmed_curve (trim validation) |
| Geom2dAPI_InterCurveCurve_Test.cxx | 3 | 119 | done (2 not ported) |
| Geom2d_AxisPlacement_Test.cxx | 14 | 185 | not ported (handle wrapper of gp_Ax2d) |
| Geom2d_BezierCurve_Test.cxx | 35 | 495 | done (copy/WeightsArray tests not ported) |
| Geom2d_BSplineCurve_Test.cxx | 32 | 505 | done (copy/WeightsArray tests not ported) |
| Geom2d_CartesianPoint_Test.cxx | 13 | 136 | not ported (handle wrapper of gp_Pnt2d) |
| Geom2d_Circle_Test.cxx | 18 | 177 | done (in geom2d_conics) |
| Geom2d_CurveEval_Test.cxx | 7 | 184 | done |
| Geom2d_Direction_Test.cxx | 15 | 159 | not ported (handle wrapper of gp_Dir2d) |
| Geom2d_Ellipse_Test.cxx | 17 | 173 | done (in geom2d_conics) |
| Geom2dEval_AHTBezierCurve_Test.cxx | 29 | 439 | not ported (non-core curve family) |
| Geom2dEval_ArchimedeanSpiralCurve_Test.cxx | 13 | 187 | not ported (non-core curve family) |
| Geom2dEval_CircleInvoluteCurve_Test.cxx | 14 | 193 | not ported (non-core curve family) |
| Geom2dEval_LogarithmicSpiralCurve_Test.cxx | 13 | 207 | not ported (non-core curve family) |
| Geom2dEval_SineWaveCurve_Test.cxx | 12 | 175 | not ported (non-core curve family) |
| Geom2dEval_TBezierCurve_Test.cxx | 27 | 478 | not ported (non-core curve family) |
| Geom2dGcc_Circ2d2TanOn_Test.cxx | 1 | 75 | done (in geom2d_gcc) |
| Geom2dGcc_Circ2d2TanRad_Test.cxx | 1 | 81 | done (in geom2d_gcc) |
| Geom2dGridEval_BezierCurve_Test.cxx | 6 | 190 | not ported (batch-evaluation internals) |
| Geom2dGridEval_Curve_Test.cxx | 29 | 680 | not ported (batch-evaluation internals) |
| Geom2dGridEval_Ellipse_Test.cxx | 5 | 157 | not ported (batch-evaluation internals) |
| Geom2dGridEval_Hyperbola_Test.cxx | 5 | 142 | not ported (batch-evaluation internals) |
| Geom2dGridEval_Parabola_Test.cxx | 6 | 160 | not ported (batch-evaluation internals) |
| Geom2dHash_CurveHasher_Test.cxx | 39 | 654 | not ported (hashing internals) |
| Geom2d_Hyperbola_Test.cxx | 14 | 152 | done (in geom2d_conics) |
| Geom2d_Line_Test.cxx | 20 | 188 | done |
| Geom2d_OffsetCurve_Test.cxx | 5 | 129 | done |
| Geom2d_Parabola_Test.cxx | 13 | 141 | done (in geom2d_conics) |
| Geom2d_Transformation_Test.cxx | 19 | 284 | done (Form() checks not ported) |
| Geom2d_TrimmedCurve_Test.cxx | 17 | 197 | done |
| Geom2d_VectorWithMagnitude_Test.cxx | 24 | 222 | not ported (handle wrapper of gp_Vec2d) |

## TKG3d (45 files)

| file | tests | lines | status |
|---|---|---|---|
| GeomAdaptor_Curve_Test.cxx | 17 | 288 | merged into geom_elementary (trim validation) |
| GeomAdaptor_TransformedCurve_Test.cxx | 16 | 366 | not ported (adaptor internals) |
| GeomAdaptor_TransformedSurface_Test.cxx | 12 | 300 | not ported (adaptor caching internals) |
| GeomAPI_ExtremaCurveCurve_Test.cxx | 47 | 933 | done (sorted-by-distance; see header) |
| GeomAPI_Interpolate_Test.cxx | 1 | 72 | done |
| Geom_BezierCurve_Test.cxx | 38 | 557 | done (copy/WeightsArray not ported) |
| Geom_BezierSurface_Test.cxx | 45 | 729 | done (copy/WeightsArray not ported) |
| Geom_BSplineCurve_Test.cxx | 43 | 736 | done (copy/WeightsArray not ported) |
| Geom_BSplineSurface_Test.cxx | 50 | 1053 | done (copy/WeightsArray not ported) |
| Geom_Circle_Test.cxx | 8 | 135 | done (in geom_elementary) |
| Geom_CurveEval_Test.cxx | 11 | 257 | done |
| GeomEval_AHTBezierCurve_Test.cxx | 18 | 344 | not ported (non-core family) |
| GeomEval_AHTBezierSurface_Test.cxx | 14 | 314 | not ported (non-core family) |
| GeomEval_CircularHelicoidSurface_Test.cxx | 14 | 210 | not ported (non-core family) |
| GeomEval_CircularHelixCurve_Test.cxx | 14 | 235 | not ported (non-core family) |
| GeomEval_EllipsoidSurface_Test.cxx | 20 | 459 | not ported (non-core family) |
| GeomEval_HyperboloidSurface_Test.cxx | 19 | 346 | not ported (non-core family) |
| GeomEval_HypParaboloidSurface_Test.cxx | 19 | 357 | not ported (non-core family) |
| GeomEval_ParaboloidSurface_Test.cxx | 16 | 285 | not ported (non-core family) |
| GeomEval_SineWaveCurve_Test.cxx | 11 | 169 | not ported (non-core family) |
| GeomEval_TBezierCurve_Test.cxx | 18 | 336 | not ported (non-core family) |
| GeomEval_TBezierSurface_Test.cxx | 17 | 388 | not ported (non-core family) |
| GeomGridEval_BezierCurve_Test.cxx | 10 | 295 | not ported (batch evaluation) |
| GeomGridEval_BezierSurface_Test.cxx | 14 | 500 | not ported (batch evaluation) |
| GeomGridEval_BSplineSurface_Test.cxx | 18 | 622 | not ported (batch evaluation) |
| GeomGridEval_Cone_Test.cxx | 5 | 186 | not ported (batch evaluation) |
| GeomGridEval_Curve_Test.cxx | 39 | 1026 | not ported (batch evaluation) |
| GeomGridEval_Cylinder_Test.cxx | 5 | 186 | not ported (batch evaluation) |
| GeomGridEval_Ellipse_Test.cxx | 4 | 141 | not ported (batch evaluation) |
| GeomGridEval_Hyperbola_Test.cxx | 4 | 124 | not ported (batch evaluation) |
| GeomGridEval_OffsetSurface_Test.cxx | 14 | 447 | not ported (batch evaluation) |
| GeomGridEval_Parabola_Test.cxx | 4 | 130 | not ported (batch evaluation) |
| GeomGridEval_Sphere_Test.cxx | 5 | 186 | not ported (batch evaluation) |
| GeomGridEval_SurfaceOfExtrusion_Test.cxx | 8 | 319 | not ported (batch evaluation) |
| GeomGridEval_SurfaceOfRevolution_Test.cxx | 7 | 258 | not ported (batch evaluation) |
| GeomGridEval_Surface_Test.cxx | 27 | 880 | not ported (batch evaluation) |
| GeomGridEval_Torus_Test.cxx | 5 | 186 | not ported (batch evaluation) |
| GeomHash_CurveHasher_Test.cxx | 40 | 692 | not ported (hashing internals) |
| GeomHash_MeshHasher_Test.cxx | 10 | 194 | not ported (hashing internals) |
| GeomHash_SurfaceHasher_Test.cxx | 43 | 754 | not ported (hashing internals) |
| Geom_Line_Test.cxx | 10 | 119 | done (in geom_elementary) |
| Geom_OffsetCurve_Test.cxx | 4 | 106 | not ported (copy tests only) |
| Geom_OffsetSurface_Test.cxx | 33 | 1258 | done (copy tests not ported; three-way merged) |
| Geom_Plane_Test.cxx | 8 | 111 | done (in geom_elementary) |
| Geom_SurfaceEval_Test.cxx | 9 | 280 | done |

## TKGeomBase (69 files)

| file | tests | lines | status |
|---|---|---|---|
| AdvApp2Var_Context_Test.cxx | 1 | 66 | not ported |
| AdvApp2Var_Framework_Test.cxx | 3 | 127 | not ported |
| AdvApp2Var_Iso_Test.cxx | 1 | 38 | not ported |
| AdvApp2Var_Network_Test.cxx | 2 | 114 | not ported |
| AdvApp2Var_Node_Test.cxx | 1 | 33 | not ported |
| AppCont_ContMatrices_Test.cxx | 9 | 279 | not ported |
| BndLib_Test.cxx | 63 | 1482 | done (in bndlib; duplicate adaptor cases merged) |
| Extrema_ExtPC_Test.cxx | 1 | 61 | done (in extrema_ext_pc) |
| ExtremaPC2d_BezierCurve_Test.cxx | 8 | 153 | done (in extrema_pc2d_curves) |
| ExtremaPC2d_BSplineCurve_Test.cxx | 7 | 144 | done (in extrema_pc2d_curves) |
| ExtremaPC2d_Circle_Test.cxx | 12 | 163 | done (in extrema_pc2d_conics) |
| ExtremaPC2d_Curve_Test.cxx | 10 | 231 | partial (in extrema_pc2d_general; see Not ported) |
| ExtremaPC2d_DistanceFunction_Test.cxx | 4 | 104 | not ported |
| ExtremaPC2d_Ellipse_Test.cxx | 9 | 134 | done (in extrema_pc2d_conics) |
| ExtremaPC2d_GridEvaluator_Test.cxx | 7 | 146 | not ported |
| ExtremaPC2d_Hyperbola_Test.cxx | 9 | 146 | done (in extrema_pc2d_curves) |
| ExtremaPC2d_Line_Test.cxx | 8 | 113 | done (in extrema_pc2d_conics) |
| ExtremaPC2d_OffsetCurve_Test.cxx | 7 | 140 | done (in extrema_pc2d_general; cache test not ported) |
| ExtremaPC2d_OtherCurve_Test.cxx | 5 | 131 | not ported |
| ExtremaPC2d_Parabola_Test.cxx | 10 | 141 | done (in extrema_pc2d_curves) |
| ExtremaPC2d_Test.cxx | 6 | 154 | partial (in extrema_pc2d_general; internal helpers not ported) |
| ExtremaPC_BezierCurve_Test.cxx | 31 | 772 | done (in extrema_pc_bezier) |
| ExtremaPC_BSplineCurve_Test.cxx | 27 | 793 | done (in extrema_pc_bspline; grid-utility test not ported) |
| ExtremaPC_Circle_Test.cxx | 26 | 542 | done (in extrema_pc_circle) |
| ExtremaPC_Comparison_Test.cxx | 17 | 640 | done (in extrema_pc_comparison; reference is extrema::point_curve) |
| ExtremaPC_Curve_Test.cxx | 47 | 1225 | done (in extrema_pc_curve; see Not ported) |
| ExtremaPC_Ellipse_Test.cxx | 30 | 615 | done (in extrema_pc_ellipse) |
| ExtremaPC_ExtendedGeometry_Test.cxx | 40 | 1025 | done (in extrema_pc_extended; reference is extrema::point_curve) |
| ExtremaPC_Hyperbola_Test.cxx | 24 | 490 | done (in extrema_pc_hyperbola) |
| ExtremaPC_Line_Test.cxx | 31 | 555 | done (in extrema_pc_line) |
| ExtremaPC_OffsetCurve_Test.cxx | 13 | 388 | done (in extrema_pc_offset; scaled-adaptor test not ported) |
| ExtremaPC_Parabola_Test.cxx | 24 | 488 | done (in extrema_pc_parabola) |
| ExtremaPC_SearchMode_Test.cxx | 30 | 1014 | done (in extrema_pc_search_mode; reference is extrema::point_curve) |
| gce_MakeCirc2d_Test.cxx | 3 | 60 | done (in construct) |
| gce_MakeCone_Test.cxx | 9 | 169 | done (in construct) |
| gce_MakeCylinder_Test.cxx | 3 | 69 | done (in construct) |
| gce_MakeElips_Test.cxx | 3 | 58 | done (in construct) |
| gce_MakeHypr_Test.cxx | 24 | 279 | done (in construct) |
| GC_MakeArcOfCircle_Test.cxx | 3 | 53 | done (in construct) |
| GC_MakeCircle2d_Test.cxx | 3 | 53 | done (in construct) |
| GC_MakeConicalSurface_Test.cxx | 2 | 40 | done (in construct) |
| GC_MakeParabola2d_Test.cxx | 3 | 88 | done (in construct) |
| GC_MakePlane_Test.cxx | 3 | 53 | done (in construct) |
| GC_MakeSegment2d_Test.cxx | 6 | 82 | done (in construct) |
| GCPnts_AbscissaPoint_Test.cxx | 5 | 80 | done (in gcpnts) |
| GCPnts_QuasiUniformDeflection_Test.cxx | 1 | 98 | done (in gcpnts) |
| Geom2dConvert_CompCurveToBSplineCurve_Test.cxx | 1 | 58 | done (in convert) |
| GeomBndLib_Curve2d_Test.cxx | 21 | 542 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_Curve_Test.cxx | 25 | 621 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_OffsetCurve2d_Test.cxx | 10 | 281 | not ported (comparisons with BndLib only) |
| GeomBndLib_OffsetCurve_Test.cxx | 14 | 355 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_OffsetSurface_Test.cxx | 15 | 375 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_SurfaceOfExtrusion_Test.cxx | 11 | 352 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_SurfaceOfRevolution_Test.cxx | 15 | 444 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomBndLib_Surface_Test.cxx | 24 | 599 | partial (in bndlib; comparisons with BndLib not ported) |
| GeomConvert_CompCurveToBSplineCurve_Test.cxx | 8 | 378 | done (in convert) |
| GeomConvert_Test.cxx | 3 | 78 | done (in convert) |
| GeomLib_CheckCurveOnSurface_Test.cxx | 7 | 198 | done |
| GeomLib_Test.cxx | 15 | 435 | partial (see Not ported) |
| GeomLProp_CLProps2d_Test.cxx | 18 | 230 | done (in lprop) |
| GeomLProp_CurAndInf2d_Test.cxx | 13 | 221 | done (in lprop) |
| GeomProjLib_Test.cxx | 2 | 95 | done (in project) |
| GProp_PEquation_Test.cxx | 4 | 85 | done (in gprop) |
| GProp_PGProps_Test.cxx | 7 | 109 | done (in gprop) |
| Hermit_Test.cxx | 11 | 252 | not ported |
| IntAna_IntQuadQuad_Test.cxx | 10 | 323 | done (in intana) |
| LProp_CurAndInf_Test.cxx | 10 | 130 | not ported |
| ProjLib_ComputeApproxOnPolarSurface_Test.cxx | 1 | 82 | done (in project) |
| ProjLib_Cone_Test.cxx | 3 | 61 | done (in project) |

## TKBRep (61 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepAdaptor_CompCurve_Test.cxx | 4 | 252 | done |
| BRepGraph_Assembly_Test.cxx | 45 | 1420 | done (22 not ported) |
| BRepGraph_BatchOps_Test.cxx | 17 | 580 | done |
| BRepGraph_Build_Test.cxx | 59 | 1327 | done (11 not ported) |
| BRepGraph_CacheMesh_Test.cxx | 22 | 618 | not ported |
| BRepGraph_CacheRegistry_Test.cxx | 14 | 472 | not ported |
| BRepGraph_ChildExplorer_Test.cxx | 37 | 1110 | done (3 not ported) |
| BRepGraph_Compact_Test.cxx | 26 | 1287 | done (6 not ported) |
| BRepGraph_Convenience_Test.cxx | 17 | 294 | done (4 not ported) |
| BRepGraph_Copy_Test.cxx | 23 | 1067 | done (7 not ported) |
| BRepGraph_Deduplicate_Test.cxx | 82 | 2729 | done (15 not ported, repetitive cases merged) |
| BRepGraph_DeferredInvalidation_Test.cxx | 16 | 413 | not ported |
| BRepGraph_DefsIterator_Test.cxx | 11 | 267 | done (2 not ported) |
| BRepGraph_EdgeCases_Test.cxx | 12 | 283 | done (8 not ported) |
| BRepGraph_EventBus_Test.cxx | 26 | 764 | not ported |
| BRepGraph_Fuzz_Test.cxx | 2 | 224 | not ported (no tests) |
| BRepGraph_Geometry_Test.cxx | 44 | 1335 | done (2 not ported) |
| BRepGraphInc_PopulateRegression_Test.cxx | 9 | 501 | not ported |
| BRepGraphInc_Test.cxx | 80 | 2790 | partial (11 cases; import/reconstruct/relations-index tests not ported) |
| BRepGraph_ItemId_Test.cxx | 4 | 64 | not ported |
| BRepGraph_Iterator_Test.cxx | 14 | 188 | done (1 not ported) |
| BRepGraph_LayerHistory_Test.cxx | 39 | 1106 | done (6 not ported) |
| BRepGraph_LayerIterator_Test.cxx | 4 | 121 | not ported |
| BRepGraph_LayerTopoSupplement_Test.cxx | 22 | 1010 | not ported |
| BRepGraph_Lock_Test.cxx | 19 | 547 | not ported |
| BRepGraph_MutationGen_Test.cxx | 14 | 398 | done (3 not ported) |
| BRepGraph_MutGuard_Test.cxx | 11 | 251 | not ported |
| BRepGraph_NodeId_Test.cxx | 15 | 221 | done (5 not ported) |
| BRepGraph_ParentExplorer_Test.cxx | 23 | 846 | done (5 not ported) |
| BRepGraph_PermissionUpdate_Test.cxx | 15 | 538 | not ported |
| BRepGraph_Polygon_Test.cxx | 7 | 397 | done (4 not ported) |
| BRepGraph_Reconstruct_Test.cxx | 24 | 649 | not ported |
| BRepGraph_RefId_Test.cxx | 19 | 569 | done (ids + use stable ids; generation/stamp parts not ported) |
| BRepGraph_RefsIterator_Test.cxx | 8 | 229 | done (2 not ported) |
| BRepGraph_RelatedIterator_Test.cxx | 15 | 357 | done (1 not ported) |
| BRepGraph_RepId_Test.cxx | 10 | 106 | done (5 not ported) |
| BRepGraph_ReplaceVertex_Test.cxx | 6 | 212 | done |
| BRepGraph_ReverseIterator_Test.cxx | 20 | 449 | done (1 not ported) |
| BRepGraph_Reverse_Test.cxx | 16 | 466 | done (5 not ported) |
| BRepGraph_ScenarioMatrix_Test.cxx | 18 | 1399 | partial (7 cases; TopoDS reconstruct/re-import scenarios not ported) |
| BRepGraph_SeamRedesign_Test.cxx | 9 | 426 | done (2 not ported) |
| BRepGraph_ShapesViewImport_Test.cxx | 44 | 1132 | done (programmatic construction; 11 not ported) |
| BRepGraph_Sharing_Test.cxx | 17 | 396 | done |
| BRepGraph_SparseModel_Test.cxx | 28 | 754 | done (empty-container cases merged) |
| BRepGraph_Test.cxx | 95 | 2116 | done (duplicates of other modules merged; ~40 not ported) |
| BRepGraph_Tool_Test.cxx | 16 | 242 | done (3 not ported) |
| BRepGraph_Transform_Test.cxx | 19 | 913 | done (3 not ported) |
| BRepGraph_TypedIdDispatch_Test.cxx | 2 | 49 | not ported |
| BRepGraph_Validate_Test.cxx | 23 | 1180 | done (5 not ported) |
| BRepGraph_VersionStamp_Test.cxx | 22 | 302 | done (8 not ported) |
| BRepGraph_Views_Test.cxx | 81 | 1513 | partial (15 cases; view-parity, cache, UID-table, mutation-guard and TopoDS tests not ported) |
| BRepTools_GTrsfModification_Test.cxx | 13 | 491 | done (2 not ported) |
| BRepTools_ReShape_Test.cxx | 16 | 406 | done |
| BRepTools_Test.cxx | 11 | 675 | done (binary format parts dropped) |
| BRepTools_TrsfModification_Test.cxx | 4 | 177 | done (1 not ported) |
| BRep_Tool_Test.cxx | 13 | 263 | done |
| TopExp_Test.cxx | 9 | 203 | done |
| TopoDS_Builder_Test.cxx | 19 | 547 | done (3 not ported) |
| TopoDS_Edge_Test.cxx | 1 | 37 | not ported |
| TopoDS_Iterator_Test.cxx | 17 | 471 | done (2 not ported) |
| TopoDS_TShape_Test.cxx | 8 | 313 | done (3 not ported) |

## TKTopAlgo (16 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepBndLib_Test.cxx | 6 | 234 |  |
| BRepBuilderAPI_Copy_Test.cxx | 6 | 167 |  |
| BRepBuilderAPI_MakeEdge_Test.cxx | 6 | 132 |  |
| BRepBuilderAPI_MakeFace_Test.cxx | 4 | 116 |  |
| BRepBuilderAPI_MakeShell_Test.cxx | 1 | 50 |  |
| BRepBuilderAPI_MakeWire_Test.cxx | 1 | 72 |  |
| BRepBuilderAPI_Transform_Test.cxx | 5 | 135 |  |
| BRepCheck_Face_Test.cxx | 4 | 199 |  |
| BRepClass3d_SolidClassifier_Test.cxx | 5 | 83 |  |
| BRepClass_FaceExplorer_Test.cxx | 8 | 386 |  |
| BRepExtrema_DistShapeShape_Test.cxx | 3 | 94 |  |
| BRepExtrema_ShapeProximity_Test.cxx | 20 | 568 |  |
| BRepGProp_Test.cxx | 17 | 486 |  |
| BRepLib_MakeWire_Test.cxx | 1 | 25 |  |
| BRepOffsetAPI_ThruSections_Test.cxx | 3 | 393 |  |
| BRepTopAdaptor_FClass2d_Test.cxx | 5 | 202 |  |

## TKPrim (8 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepPrimAPI_MakeBox_Test.cxx | 6 | 124 | done |
| BRepPrimAPI_MakeCone_Test.cxx | 6 | 113 | done |
| BRepPrimAPI_MakeCylinder_Test.cxx | 5 | 99 | done |
| BRepPrimAPI_MakePrism_Test.cxx | 1 | 42 | done |
| BRepPrimAPI_MakeRevol_Test.cxx | 3 | 120 | done |
| BRepPrimAPI_MakeSphere_Test.cxx | 5 | 102 | done |
| BRepPrimAPI_MakeTorus_Test.cxx | 15 | 386 | done |
| BRepPrimAPI_MakeWedge_Test.cxx | 5 | 94 | done |

## TKBO (14 files)

| file | tests | lines | status |
|---|---|---|---|
| BOPAlgo_BOP_Test.cxx | 22 | 499 |  |
| BOPAlgo_CellsBuilder_Test.cxx | 2 | 185 |  |
| BOPAlgo_MakePeriodic_Test.cxx | 9 | 998 |  |
| BOPAlgo_OpenSolid_Test.cxx | 9 | 473 |  |
| BOPAlgo_PaveFiller_Test.cxx | 4 | 352 |  |
| BRepAlgoAPI_BuilderAlgo_Test.cxx | 10 | 90 |  |
| BRepAlgoAPI_Common_Test.cxx | 2 | 58 | done |
| BRepAlgoAPI_Cut_Test_1.cxx | 29 | 1409 |  |
| BRepAlgoAPI_Cut_Test.cxx | 84 | 1910 | done |
| BRepAlgoAPI_Fuse_Test.cxx | 100 | 2051 | done |
| BRepAlgoAPI_Simplify_Test.cxx | 2 | 135 |  |
| BRepAlgoAPI_Splitter_Test.cxx | 10 | 354 |  |
| BRepTools_History_Test.cxx | 6 | 359 |  |
| IntTools_FaceFace_Test.cxx | 10 | 536 |  |

## TKBool (5 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepAlgoAPI_Cut_Test.cxx | 4 | 222 |  |
| BRepAlgoAPI_Fuse_Test.cxx | 6 | 204 |  |
| BRepAlgoAPI_Section_Test.cxx | 4 | 168 |  |
| BRepFill_PipeShell_Test.cxx | 2 | 85 |  |
| TopOpeBRepDS_BuildTool_Test.cxx | 1 | 43 |  |

## TKGeomAlgo (37 files)

| file | tests | lines | status |
|---|---|---|---|
| GccAna_Circ2d3Tan_Test.cxx | 18 | 343 |  |
| Geom2dAPI_ExtremaCurveCurve_Test.cxx | 1 | 49 |  |
| Geom2dAPI_InterCurveCurve_Test.cxx | 6 | 235 | done (2 not ported) |
| Geom2dAPI_Interpolate_Test.cxx | 1 | 80 |  |
| Geom2dAPI_PointsToBSpline_Test.cxx | 2 | 56 |  |
| Geom2dConvert_BSplineCurveToBezierCurve_Test.cxx | 1 | 57 |  |
| Geom2dGcc_Circ2d2TanRad_Test.cxx | 7 | 264 | done (in geom2d_gcc) |
| Geom2dGcc_Circ2d3Tan_Test.cxx | 12 | 480 |  |
| Geom2dGcc_Lin2d2Tan_Test.cxx | 2 | 82 |  |
| Geom2dHatch_Elements_Test.cxx | 3 | 84 |  |
| Geom2dHatch_Intersector_Test.cxx | 3 | 106 |  |
| GeomAPI_ExtremaCurveSurface_Test.cxx | 6 | 172 |  |
| GeomAPI_ExtremaSurfaceSurface_Test.cxx | 2 | 48 |  |
| GeomAPI_IntCS_Test.cxx | 5 | 148 |  |
| GeomAPI_IntSS_Test.cxx | 6 | 3752 |  |
| GeomAPI_PointsToBSplineSurface_Test.cxx | 1 | 46 |  |
| GeomAPI_PointsToBSpline_Test.cxx | 1 | 40 |  |
| GeomAPI_ProjectPointOnSurf_Test.cxx | 2 | 103 |  |
| GeomFill_BSplineCurves_Test.cxx | 5 | 271 |  |
| GeomFill_CorrectedFrenet_Test.cxx | 4 | 170 |  |
| GeomFill_Gordon_Test.cxx | 64 | 2547 |  |
| GeomFill_GuideTrihedronAC_Test.cxx | 4 | 148 |  |
| GeomFill_NSections_Test.cxx | 1 | 59 |  |
| GeomPlate_BuildPlateSurface_Test.cxx | 6 | 239 |  |
| IntCurveSurface_IntersectionPoint_Test.cxx | 3 | 58 |  |
| IntCurveSurface_InterUtils_Test.cxx | 1 | 94 |  |
| IntCurveSurface_ThePolygonOfHInter_Test.cxx | 1 | 34 |  |
| IntCurveSurface_ThePolyhedronOfHInter_Test.cxx | 6 | 136 |  |
| Intf_Tool_Test.cxx | 5 | 134 |  |
| IntPatch_PolyhedronBVH_Test.cxx | 8 | 239 |  |
| IntPatch_Polyhedron_Test.cxx | 5 | 111 |  |
| IntPolyh_Intersection_Test.cxx | 3 | 124 |  |
| IntPolyh_Point_Test.cxx | 10 | 116 |  |
| IntSurf_LineOn2S_Test.cxx | 3 | 76 |  |
| IntSurf_Quadric_Test.cxx | 1 | 43 |  |
| Plate_Plate_Test.cxx | 2 | 50 |  |
| TopTrans_SurfaceTransition_Test.cxx | 2 | 66 |  |

## TKFillet (3 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepFilletAPI_MakeChamfer_Test.cxx | 6 | 299 | done (edges chosen by position, volumes checked too) |
| BRepFilletAPI_MakeFillet_Test.cxx | 13 | 768 |  |
| ChFi3d_Hatching_Test.cxx | 2 | 99 |  |

## TKOffset (4 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepBuilderAPI_Sewing_Test.cxx | 5 | 117 |  |
| BRepOffsetAPI_MakePipeShell_Test.cxx | 2 | 322 |  |
| BRepOffsetAPI_MakeThickSolid_Test.cxx | 4 | 149 |  |
| BRepOffset_MakeOffset_Test.cxx | 24 | 1099 |  |

## TKShHealing (12 files)

| file | tests | lines | status |
|---|---|---|---|
| ShapeAnalysis_CanonicalRecognition_Test.cxx | 21 | 814 |  |
| ShapeAnalysis_Edge_Test.cxx | 5 | 91 |  |
| ShapeAnalysis_FreeBounds_Test.cxx | 11 | 332 |  |
| ShapeBuild_ReShape_Test.cxx | 4 | 176 |  |
| ShapeConstruct_ProjectCurveOnSurface_Test.cxx | 29 | 999 |  |
| ShapeFix_EdgeConnect_Test.cxx | 18 | 394 |  |
| ShapeFix_FaceConnect_Test.cxx | 16 | 341 |  |
| ShapeFix_Face_Test.cxx | 4 | 158 |  |
| ShapeFix_Shape_Test.cxx | 7 | 192 |  |
| ShapeUpgrade_FaceDivide_Test.cxx | 2 | 66 |  |
| ShapeUpgrade_ShellSewing_Test.cxx | 17 | 322 |  |
| ShapeUpgrade_UnifySameDomain_Test.cxx | 9 | 739 |  |

## TKHelix (6 files)

| file | tests | lines | status |
|---|---|---|---|
| HelixBRep_BuilderHelix_Integration_Test.cxx | 10 | 418 |  |
| HelixBRep_BuilderHelix_Test.cxx | 14 | 436 |  |
| HelixGeom_BuilderHelixCoil_Test.cxx | 4 | 132 |  |
| HelixGeom_BuilderHelix_Test.cxx | 4 | 103 |  |
| HelixGeom_HelixCurve_Test.cxx | 9 | 298 |  |
| HelixGeom_Tools_Test.cxx | 3 | 134 |  |

## TKMesh (6 files)

| file | tests | lines | status |
|---|---|---|---|
| BRepMesh_BaseMeshAlgo_Test.cxx | 3 | 277 |  |
| BRepMesh_CircleTool_Test.cxx | 1 | 84 |  |
| BRepMesh_Delaun_Test.cxx | 6 | 237 |  |
| BRepMesh_DiscretAlgoFactory_Test.cxx | 9 | 145 |  |
| BRepMesh_GeomTool_Test.cxx | 1 | 86 |  |
| BRepMesh_IncrementalMesh_Test.cxx | 9 | 463 |  |


## BRL-CAD libbg (`src/libbg/tests`)

Hand translations in `tests/brlcad_libbg/`, one module per BRL-CAD test source, each test
registered in its `CMakeLists.txt` with a counterpart of the same name (`-` in names becomes
`_minus`). BRL-CAD's return codes become geometric assertions on kentron's results (closest
points, distances, which end of a segment); its default distance tolerance (0.0005) is kept
where the cases depend on it.

| source | tests | kentron |
|---|---|---|
| plane_pt.c | 48 | `Plane::through_points`, `Line::distance_to_point`, `math::closest_point_on_segment`, `Plane` parameters |
| plane_dist.c | 11 | `math::closest_point_on_segment` |
| plane_isect.c | 16 | `math::closest_points_on_segments` |
| lseg_pt.c | 6 | `math::closest_point_on_segment` |
| lseg_lseg.c | 1 | `math::closest_points_on_segments` (BRL-CAD only prints the result: checked against a dense search) |
| tri_tri_isect.c | 5 (+1) | `math::triangles_intersect` |
| tri_tri_isect_coplanar.cpp | 4 | `math::triangles_intersect` |
| tri_closest_pt.c | 2 | `math::closest_point_on_triangle` |
| tri_area.c | 1 | `math::triangle_area` |
| trimesh_area.c | 1 | `mesh::TriangleMesh::area` |

Deviation: the three `bg_plane_dist_pt3_line3_*` cases run BRL-CAD's segment test
(`bg_dist_pnt3_lseg3`) with return codes of its line test; BRL-CAD's driver accepts them only
because it fails a case when every value differs (`&&` for `||`). Their closest points and
distances are checked; their return codes are not.

Not ported (no kentron counterpart): the `bg_tri_tri_isect_coplanar_*_area` cases (an overlap
of non-zero area), `polygon_triangulate.c`, `poly2tri.c`, `polygon_op.c`, `chull.c`,
`chull3d.c`, `ballpivot.c`, `obr.c`, `sat.c`, `tri_ray_isect.cpp`, and the triangle mesh
repairs (`trimesh_sync.c`, `trimesh_repair.cpp`, `trimesh_remesh.cpp`).
