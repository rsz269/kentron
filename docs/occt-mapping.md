# OCCT → kentron mapping

kentron's functionality is modelled on OCCT, but its names and structure are its own. This
table maps OCCT concepts to kentron ones; it grows as tests are ported. Types also carry
`#[doc(alias = "...")]` entries, so searching the rustdoc for an OCCT name finds them.

## Types

| OCCT | kentron | notes |
|---|---|---|
| `gp_Pnt`, `gp_Pnt2d` | `math::Point`, `Point2` | public `x`, `y`, `z` fields |
| `gp_Vec`, `gp_Vec2d` | `math::Vector`, `Vector2` | |
| `gp_XYZ`, `gp_XY` | `math::Vector`, `Vector2` | no separate coordinate-triple type |
| `gp_Dir`, `gp_Dir2d` | `math::Direction`, `Direction2` | always unit length; `Direction::X`, `NEG_Z`, ... replace `gp_Dir::D` and `gp::DX()` |
| `gp_Ax1`, `gp_Ax2d` | `math::Axis`, `Axis2` | `Axis::X/Y/Z` |
| `gp_Ax2`, `gp_Ax3` | `math::Frame` | one type; handedness is a property (`is_right_handed`) |
| `gp_Mat` | `math::Matrix3` | indexed `m[(row, col)]` from 0 |
| `gp_Trsf`, `gp_Trsf2d` | `math::Transform`, `Transform2` | built by constructors (`Transform::rotation(...)`); `a * b` applies `b` first |
| `gp_GTrsf` | `math::AffineTransform` | |
| `gp_Quaternion`, `gp_EulerSequence` | `math::Quaternion`, `EulerSequence` | |
| `gp_Lin` + `Geom_Line` | `geometry::curve::Line` | gp and Geom variants merged |
| `gp_Lin2d` + `Geom2d_Line` | `geometry::curve2::Line2` | |
| `gp_Circ` + `Geom_Circle` | `geometry::curve::Circle` | |
| `gp_Pln` + `Geom_Plane` | `geometry::surface::Plane` | |
| `gp_Torus` + `Geom_ToroidalSurface` | `geometry::surface::Torus` | |
| `Geom_Curve`, `Geom2d_Curve`, `Geom_Surface` | `geometry::curve::Curve`, `curve2::Curve2`, `surface::Surface` | closed enums; `DownCast` is a `match` |
| `Handle(T)` | `T` (owned value) or a `topology::Id<T>` | no reference counting |
| `Precision::Confusion()`, `Angular()`, `PConfusion()`, `Infinite()` | `tolerance::LINEAR`, `ANGULAR`, `PARAMETRIC`, `INFINITE` | constants |
| `Poly_Triangulation`, `Poly_Polygon3D`, `Poly_Polygon2D`, `Poly_PolygonOnTriangulation` | `mesh::TriangleMesh`, `Polyline`, `Polyline2`, `MeshPolyline` | node indices are 0-based |
| `Geom_CylindricalSurface`, `Geom_ConicalSurface`, `Geom_SphericalSurface`, `Geom_RectangularTrimmedSurface` | `geometry::surface::Cylinder`, `Cone`, `Sphere`, `TrimmedSurface` | |
| `Geom_BezierCurve`, `Geom_BezierSurface`, `Geom_TrimmedCurve` | `geometry::curve::BezierCurve`, `surface::BezierSurface`, `curve::TrimmedCurve` | |
| `BRepGraph` | `topology::Model` | the arena owning an assembly |
| `TopoDS_Shape` | `topology::Shape` | node + orientation + placement, `Copy` |
| `TopoDS_TVertex` ... `TopoDS_TCompound` | `topology::Vertex` ... `Compound` | the BRepGraph "defs" |
| BRepGraph coedge | `topology::CoEdge` | |
| BRepGraph refs (`FaceRef`, `ChildRef`, ...) | `topology::FaceUse`, `CompoundChild`, ... | |
| BRepGraph reps (`EdgeCurve3DRep`, `CoEdgeCurve2DRep`, `FaceSurfaceRep`, ...) | `topology::EdgeCurve`, `PCurve`, `FaceSurface`, ... | |
| BRepGraph product / occurrence | `topology::Product`, `Instance` | |
| `BRepGraph_UID` | `topology::StableId` | |
| `TopAbs_ShapeEnum`, `TopAbs_Orientation` | `topology::NodeKind`, `Orientation` | |
| `TopLoc_Location` | `math::Transform` | placements are plain transforms |

## Geometry

| OCCT | kentron |
|---|---|
| `Geom_Curve` / `Geom2d_Curve` (D0..DN, `FirstParameter`, `IsClosed`, `Continuity`, `ReversedParameter`, `Resolution`) | `geometry::ParametricCurve` (`point`, `d1`..`d3`, `dn`, `first_parameter`, `is_closed`, `continuity`, `reversed_parameter`, `resolution`) |
| `Geom_Curve::TransformedParameter`, `Geom2d_Curve::TransformedParameter` | `Curve::parameter_scale`, `Curve2::parameter_scale` (the factor parameters scale by) |
| `Geom_Surface` (D0..DN, `Bounds`, `UIso`/`VIso`, `UReverse`) | `geometry::ParametricSurface` (`point`, `d1`, `d2 -> SurfaceD2`, `d3 -> SurfaceD3`, `dn`, `bounds`, `u_iso`/`v_iso`, `u_reverse`) |
| `GeomAbs_Shape` | `geometry::Continuity` |
| `gp_Ax22d` | `math::Frame2` |
| `Geom_Line`/`Circle`/`Ellipse`/`Hyperbola`/`Parabola` (and `gp_*`) | `curve::Line`, `Circle`, `Ellipse`, `Hyperbola`, `Parabola` |
| `Geom2d_*` conics | `curve2::Line2`, `Circle2`, `Ellipse2`, `Hyperbola2`, `Parabola2` |
| `Geom_BezierCurve` / `Geom2d_BezierCurve` | `BezierCurve` / `BezierCurve2` (= `BezierCurve<Point2>`), 0-based poles |
| `Geom_BSplineCurve` / `Geom2d_BSplineCurve` | `BSplineCurve` / `BSplineCurve2`, 0-based poles and knots; `new_periodic`, `with_weights` |
| `Geom_BSplineSurface` | `surface::BSplineSurface` with a `KnotVector` per direction |
| `Geom_TrimmedCurve` / `Geom2d_TrimmedCurve` | `TrimmedCurve<C>` / `TrimmedCurve2` |
| `Geom_OffsetCurve` / `Geom2d_OffsetCurve` | `curve::OffsetCurve` / `curve2::OffsetCurve2` |
| `Geom_SurfaceOfLinearExtrusion`, `Geom_SurfaceOfRevolution`, `Geom_OffsetSurface` (+ `Surface()`) | `surface::ExtrusionSurface`, `RevolutionSurface`, `OffsetSurface` (+ `equivalent_surface()`) |
| `Eccentricity`, `Focal`, `Focus1/2`, `Parameter`, `Asymptote1/2`, `Directrix` | `eccentricity`, `focal_distance`, `foci`, `semi_latus_rectum`, `asymptotes`, `directrix` |
| `GeomAPI_Interpolate`, `GeomAPI_PointsToBSpline` (and 2D) | `fit::interpolate` / `interpolate_with_end_tangents` / `interpolate_with_parameters`, `fit::approximate(_with)` |
| `GeomAPI_ExtremaCurveCurve` / `Geom2dAPI_ExtremaCurveCurve` (`IsParallel`, `LowerDistance`, `NearestPoints`), `GeomAPI_ProjectPointOnCurve` | `extrema::curve_curve` / `curve_curve2` (`CurveCurveExtrema::Parallel`, `min_distance`), `extrema::closest_points`, `extrema::point_curve` / `point_curve2` |
| `GeomAPI_ExtremaCurveSurface` | `extrema::curve_surface` |
| `IsCNu` / `IsCNv`, `LocalD0`, `UKnot`/`UMultiplicity` | `is_cn_u` / `is_cn_v`, `local_point`, `u_knot`/`u_multiplicity` |
| `Geom2dAPI_InterCurveCurve` | `intersect::curves2` |
| `GCPnts_AbscissaPoint` (`Length`) | `measure::length`, `length_between`, `parameter_at_length` (and `*2`) |
| `Geom2dGcc_Circ2d2TanRad`, `Circ2dTanCen`, `Circ2d2TanOn`, `GccEnt_Position` | `tangency::with_radius`, `with_center`, `centered_on`, `tangency::Qualifier` |
| `gce_Make*` / `GC_Make*` / `GCE2d_Make*` (`gce_ErrorType`) | constructors on the types, e.g. `Circle2::through_points`, `Cone::offset`/`parallel_through`/`from_points_radii`, `Hyperbola::try_new`/`from_points`, `Circle::arc_through`/`arc_tangent`, `Line2::segment`, `Plane::offset`; errors are `construct::ConstructionError` |
| `gp::Resolution` | `tolerance::RESOLUTION` |
| `GCPnts_QuasiUniformDeflection` | `sample::quasi_uniform_deflection` |
| `GeomConvert::CurveToBSplineCurve`, `SurfaceToBSplineSurface`, `Geom2dConvert::CurveToBSplineCurve` | `convert::curve_to_bspline`, `surface_to_bspline`, `curve2_to_bspline` (trait `convert::ToBSpline`) |
| `GeomConvert_ApproxCurve` / `Geom2dConvert_ApproxCurve` | `convert::approximate_curve` / `approximate_curve2` (`CurveApproximationOptions`) |
| `GeomConvert_CompCurveToBSplineCurve` / `Geom2dConvert_CompCurveToBSplineCurve` | `convert::Concatenation<P>` (`add`, `result`, `clear`) |
| `GeomProjLib::ProjectOnPlane` | `project::curve_on_plane` |
| `ProjLib_Cone`, `ProjLib_Cylinder`, ... | `project::analytic_pcurve` |
| `ProjLib_ComputeApproxOnPolarSurface`, `ProjLib_ComputeApprox` | `project::pcurve`, `project::pcurve_on_range`, `project::polyline_pcurve_on_range` |
| `ElSLib::Parameters` | `parameters` on `Plane`, `Cylinder`, `Cone`, `Sphere`, `Torus`, and `Surface::parameters` |
| `GeomLib_CheckCurveOnSurface` (`MaxDistance`, `MaxParameter`) | `measure::curve_on_surface_deviation` |
| `GeomLib::ExtendSurfByLength`, `ExtendCurveToPoint`, `AdjustExtremity` | `extend::surface_by_length` (`SurfaceSide`), `extend::curve_to_point`, `extend::adjust_ends` |
| `GeomLib::To3d` | `Curve2::to_3d` |
| `GeomLib::IsBSplUClosed`, `IsBzUClosed` | `BSplineSurface::is_u_closed_between`, `BezierSurface::is_u_closed_between` |
| `GeomLib::NormEstim` | `props::surface_normal` (`SurfaceNormal::Regular` / `Singular`) |
| `GeomLProp_CLProps2d` | `props::curve_point2` (`CurvePoint2`; undefined tangent/normal/center are `None`) |
| `GeomLProp_CurAndInf2d` (`Perform`, `PerformInf`, `PerformCurExt`), `LProp_CIType` (`MinCur`/`MaxCur`) | `props::curvature_features2`, `inflections2`, `curvature_extrema2`; `CurvatureFeatureKind::MinRadius`/`MaxRadius` |
| `GProp_PGProps` (`Barycentre`), `GProp_PEquation` | `properties::points`, `point_masses`, `barycenter`; `properties::point_set_shape` (`PointSetShape`) |
| `Bnd_Box` (`IsVoid`, `IsOpen*`, `Enlarge`), `Bnd_Box2d` | `math::BoundingBox` (`EMPTY`, `is_empty`, infinite sides, `enlarged`), `BoundingBox2` |
| `BndLib::Add`, `BndLib_Add3dCurve`, `BndLib_Add2dCurve`, `BndLib_AddSurface`, `GeomBndLib_*` (`Add` / `AddOptimal`) | `bounds::curve`, `curve_between`, `curve2*`, `surface`, `surface_within` / `*_tight` |
| `ExtremaPC_Curve`, `ExtremaPC2d_Curve` and the per-type `ExtremaPC_*` evaluators (`Perform`, `PerformWithEndpoints`, `SearchMode`, `Domain1D`, `Result`) | `extrema::PointCurveSearch` (`new`, `with_domain`, `find`, `find_with_endpoints`, `find_with`), `SearchMode`, `PointCurveExtrema` (`min_index`, `infinite_square_distance`), `ExtremaError` |
| `IntAna_QuadQuadGeo`, `IntAna_IntQuadQuad` (`HasNextCurve`/`NextCurve`), `IntAna2d_AnaIntersection` (circles) | `intersect::quadrics` (`QuadricIntersection`), `intersect::cylinder_quadric` (`QuadricCurves::next`), `intersect::circles2` |

## Algorithms

| OCCT | kentron |
|---|---|
| `BRepPrimAPI_MakeBox`, `MakeWedge`, `MakeCylinder`, `MakeCone`, `MakeSphere`, `MakeTorus` | `primitives::Cuboid`, `Wedge`, `Cylinder`, `Cone`, `Sphere`, `Torus` (`::new(..).with_*(..).build(&mut model)`) |
| `BRepBuilderAPI_MakeVertex`, `MakeEdge`, `MakeWire` | `make::vertex`, `make::edge` / `edge_between` / `edge_from_curve`, `make::wire` |
| `BRepPrimAPI_MakePrism`, `MakeRevol` | `sweep::extrude`, `sweep::revolve` (returning a `Sweep`) |
| DRAW `profile`, `BOPTest_Utilities::CreateProfile` | `sketch::profile` with `sketch::Step`s: `F` `Start`, `O` `Origin`, `P` `Plane`, `S` `Frame` (the face's plane), `X`/`Y` `AlongX`/`AlongY`, `XX`/`YY` `ToX`/`ToY`, `L` `Forward`, `IX`/`IY` `ForwardToX`/`ForwardToY`, `T` `By`, `TT` `To`, `R` `Turn`, `RR` `Heading`, `D` `Direction`, `C` `Arc`; `W`/`WW` are `Outcome::ClosedWire`/`OpenWire` |
| `BRepOffsetAPI_MakePipe` (corrected Frenet trihedron, `BRepFill_SectionPlacement`, the `Modified` transition at corners); DRAW `pipe` | `sweep::pipe` (straight spines extrude, circular ones revolve, others sweep B-spline surfaces along a rotation-minimizing frame; the profile stands where the spine crosses its plane; a closed spine gives no ends) |
| `BRepOffsetAPI_MakePipeShell` (`SetMode` corrected Frenet, `Add(profile, WithContact, WithCorrection)`, `SetTransitionMode`, `MakeSolid`); DRAW `mksweep`, `setsweep -CF/-DT/-FR`, `addsweep [-T] [-R]`, `buildsweep [-M\|-C] [-S]` | `sweep::pipe_with(model, profile, spine, &PipeOptions { contact, correction })`; right corners (`BRepBuilderAPI_RightCorner`): `sweep::pipe_cornered(model, profile, spine, options, Corner::Mitred)` (spines of lines and arcs, closed planar profiles); several sections: `sweep::pipe_through(model, sections, spine, solid)` (sections blended linearly along the spine; guides, fixed and constant-binormal modes and round corners not supported) |
| `BRepBuilderAPI_MakeShape::Generated` | `Sweep::generated` |
| `BRepGProp::LinearProperties`, `SurfaceProperties`, `VolumeProperties` | `properties::length`, `area`, `volume` |
| `GProp_GProps::Mass`, `CentreOfMass` | `MassProperties::mass`, `center_of_mass` |
| `BRepCheck_Analyzer::IsValid` | `validity::check(&model, shape).is_valid()` |
| `BOPAlgo_CheckerSI` (DRAW `bopcheck`; `BOPAlgo_ArgumentAnalyzer`, `bopargcheck`, in part) | `validity::self_interferences(&mut model, shape)` -> `Vec<Interference { first, second }>` (vertices, edges and faces meeting without sharing a vertex; solids not checked) |
| `TopExp_Explorer` | `Model::explore` |
| `TopExp::MapShapes` | `Model::unique_subshapes` |
| `BRep_Tool::Surface(face)` | `Model::face_surface(face_id)` |
| `TopoDS::Face(shape)` (downcast) | `FaceId::try_from(shape.node)` |
| `IsDone()` / `IsNull()` checks | `Result` / (a `Shape` is never null) |
| `TopoDS_Builder::MakeWire/MakeCompound/...`, `BRep_Builder::MakeVertex/MakeEdge/MakeFace` (empty) | `Model::new_empty(kind)` |
| `TopoDS_Builder::Add` / `Remove` | `Model::add` / `Model::remove` |
| `TopoDS_TShape::EmptyCopy` | `Model::empty_copy` |
| `BRep_Builder::UpdateVertex/UpdateEdge/UpdateFace/Range` | `Model::set_vertex_point`, `set_edge_curve`, `set_edge_polyline`, `set_face_surface`, `set_face_mesh`, `set_edge_range` |
| `TopoDS_Iterator` (cumulative / not) | `Model::children` / `Model::children_local` |
| `TopExp::FirstVertex/LastVertex/Vertices`, `CommonVertex`, `MapShapesAndAncestors` | `Model::edge_vertices`, `common_vertex`, `subshape_ancestors` |
| `BRep_Tool::Pnt`, `Tolerance`, `Curve`, `Range`, `IsClosed`, `Degenerated`, `CurveOnPlane`, `Triangulation`, `Polygon3D` | `Model::point`, `tolerance`, `edge_curve`, `edge_range`, `is_closed_edge`, `is_degenerate`, `pcurve_on_plane`, `face_mesh`, `edge_polyline` |
| `BRepGraph_Iterator<T>` (`FaceIterator`, ...), full traversal | `model.table::<T>().iter()`, `.iter_all()` |
| `BRepGraph_DefsIterator` / `RefsIterator` (single-level children) | entity fields: `model[shell].faces`, then `model[face_use].face` |
| `BRepGraph_ReverseIterator`, `BRepGraph_ParentExplorer` | `Model::parents(node)`, `Model::ancestors_of::<T>(node)` |
| `BRepGraph_ChildExplorer` (avoid kind, emit boundary, cumulative loc/ori, `NodeOf`) | `Model::explore` / `explore_all` + `.skip_below`, `.stop_at`, `.compose_placement`, `.compose_orientation`, `Explorer::ancestor` |
| `BRepGraph_ParentExplorer` (`LeafLocation`, `LeafOrientation`) | `Model::explore_up` + `.of_kind`, `.direct_only`, `.skip_above`, `.stop_at`, `leaf_placement`, `leaf_orientation` |
| `EditorView` `Edges().Add`, `CoEdges().Add`, `Wires().Add` | `Model::add_edge`, `add_coedge`, `add_wire` |
| `Wires().Reverse`, `SetCoEdgeOrder`, `RemoveCoEdge`, append / replace edge (and their `Check*` dry runs) | `Model::reverse_wire`, `set_coedge_order`, `remove_coedge`, `append_coedge`, `replace_wire_edge` (dry run: try on `model.clone()`) |
| `Edges().Reverse`, `Edges().ReplaceVertex` | `Model::reverse_edge`, `replace_edge_vertex` |
| `CoEdges().SetPCurve`, `SetParamRange`, `ClearPCurve`, `SetOrientation`; `Faces().ClearSurface`, `Edges().ClearCurve` | `Model::set_pcurve(_with_range)`, `set_coedge_range`, `clear_pcurve`, `set_coedge_orientation`, `clear_face_surface`, `clear_edge_curve` |
| `BRepGraph_Tool::CoEdge::PCurve`, `Range`, `PCurveAdaptor`, `SameRange`, `SameParameter`, `SeamPair`; `Edge::IsSeamOnFace` | `Model::pcurve`, `coedge_range`, `effective_pcurve`, `is_same_range`, `is_same_parameter`, `seam_pair`, `is_seam_on_face` |
| batch `Append` / `RemoveFaces` / `ReplaceChild` / `Products().Append(array)` | `Model::add_all`, `remove_uses`, `replace_child`, `add_instances`, `remove_instance` |
| `BRepGraph_Compact::Perform` (+ `Result`) | `Model::compact` -> `CompactReport` |
| `BRepGraph_Deduplicate::Perform` (+ `Options{AnalyzeOnly, HistoryMode, MergeEntitiesWhenSafe, CompTolerance, HashTolerance}`) | `Model::deduplicate` / `deduplicate_with(DeduplicateOptions { analyze_only, record_history, merge_entities, angular_tolerance, linear_tolerance })` |
| `BRepGraph_Deduplicate::Result` (`NbCanonicalSurfaces/Curves`, `NbSurface/CurveRewrites`, `NbMerged*`, `IsEntityMergeApplied`, `NbHistoryRecords`) | `DeduplicateReport` (`distinct_surfaces/curves`, `shared_surfaces/curves`, `merged(kind)`, `merge_applied`, `history_events`) |
| `BRepGraph::ShapesView::Add` of a copied shape (`BRepBuilderAPI_Copy` into a fresh graph) | `Model::import(&source, shape)` |
| `BRepGraph_LayerHistory` (`Record`, `RecordBatch`, `RecordDeleted`, `RecordReplaced`, `FindModified`, `FindGenerated`, `FindOriginals`, `FindOriginal`, `FindDerived`, `IsDeleted`); `*Uid` variants | `topology::History` (`record`/`record_as`, `record_pairs`, `record_deleted`, `record_replaced`, `modified`, `generated`, `originals`, `find_original`, `find_derived`, `is_deleted`); `History<StableId>` |
| `Editor().Gen().ApplyModification` | `Model::apply_modification` |
| `LayerRegistry().Ensure<BRepGraph_LayerHistory>()` / `Find<...>()` | `Model::enable_history()` / `Model::history()` (opt-in; `None` until enabled) |
| `Editor().Gen().ReplaceNode`, `Edges().Split`, `Edges().SetTolerance`, `Faces().SetTolerance`, `Mut(..).MarkDirty()` | `Model::replace_node`, `split_edge`, `set_edge_tolerance`, `set_face_tolerance`, `touch` |
| `UIDs().Of`, `NodeIdFrom`, `RefIdFrom` | `Model::stable_id` / `node_stable_id`, `find_stable_node`, `find_stable::<T>` |
| `OwnGen` / `SubtreeGen` | `Model::revision(node)` -> `Revision { own, subtree }` |
| `UIDs().StampOf`, `IsStale`, `VersionStamp::IsSameItem` | `Model::stamp`, `is_stale`, `Stamp::is_same_item` |
| `BRepBuilderAPI_Copy` | `modify::copy` |
| `BRepGraph_Copy::Perform` / `CopyNode` | `Model::clone` / `Model::import` |
| `BRepGraph_Transform::Perform` / `TransformNode` / `MoveRef` | `modify::transform_model` / `import` + `modify::transform_in_place` / `Model::move_compound_child`, `move_instance` |
| `BRepGraph_Tool::Wire::IsClosed` | `Model::is_closed_wire` |
| `BRepGraph_Tool::Edge::IsBoundary` / `IsManifold` | `Model::is_boundary_edge` / `is_manifold_edge` |
| `BRepGraph_Validate` `Lightweight` / `Audit` | `validity::check_model_structure` / `check_model` |
| `BRepGraph_Validate::Issue` (`Sev`, `Description`), `NbIssues(sev)` | `validity::Problem` (`severity`, `kind: ProblemKind`, `description`), `Report::errors()` / `warnings()` / `has(kind, severity)` |
| `BRepGraph_RelatedIterator` (+ `RelationKind`) | `Model::related(node)` (+ `topology::Relation`) |
| `BRepGraph_Tool::Face::OuterWire`, `Wire::FaceOf`, `Edge::StartVertexId/EndVertexId`, `Edge::FindByVertices`, `Edge::FindPCurveCoEdgeId`, `Shell::IsClosed`, `CoEdge::IsSeam` | `Model::outer_wire`, `face_of_wire`, `edge_end_vertices`, `find_edge`, `find_coedge` / `find_coedge_oriented`, `is_closed_shell`, `is_seam` |
| `BRepGraph::EditorView` `Vertices().Add`, `Gen().RemoveNode` (+ `CleanupRemovedReferences`), `Gen().RemoveRef` | `Model::add_vertex`, `remove_node`, `remove_use` |
| `Products().Add(root)` / `Add()`, `Products().Append`, `RootProductIds` | `Model::new_part`, `new_assembly`, `add_instance`, `root_products` |
| `Products().IsPart` / `IsAssembly`, `Occurrences().ParentProduct` + `OccurrenceLocation`, occurrences of a product | `Product::is_part` / `is_assembly`, `Model::instance_placement`, `Model::instances_of` |
| `Editor().Gen().RemoveSubgraph(occurrence/product)` | `Model::remove_instance` / `remove_node(product)` |
| `Editor().Gen().RemoveSubgraph(topology node)` | `Model::remove_subtree` |
| `BRepGraph_Validate` | `validity::check_model` |
| `Topo().Faces().Nb()` | `Model::count(NodeKind::Face)` |
| `BRepTools::UVBounds`, `Clean` | `Model::face_uv_bounds`, `clear_meshes` |
| `BRepTools_WireExplorer` | `Model::ordered_edges` |
| `BRepTools_ReShape` (`Replace`, `Remove`, `Value`, `ValueLeaf`, `Apply`) | `topology::Substitution` (`replace`, `remove`, `get`, `resolve`, `apply`) |
| `BRepAdaptor_CompCurve` | `topology::WireCurve` |
| `TopoDS_Shape::Reversed`, `Located`, `Moved`, `IsSame` | `Shape::reversed`, `with_placement`, `moved`, `is_same` |
| `BRepTools::Write/Read` | `io::brep::write/read` |
| `StlAPI_Writer` (DRAW `writestl`) | `io::stl::write` (nodes and triangles from `meshing::triangulate`) |
| `BRepMesh_IncrementalMesh` (DRAW `incmesh`; CadQuery `tessellate`) | `meshing::mesh`, `mesh_with(MeshOptions)`, `triangulate` |
| `BRepExtrema_DistShapeShape` (DRAW `distmini`; CadQuery `distance`) | `distance::distance` -> `ShapeDistance` |
| (a distance only needed up to a bound) | `distance::distance_within(model, a, b, limit)` |
| `BRepClass3d_SolidClassifier` (CadQuery `isInside`) | `boolean::classify_point` -> `PointLocation` |
| `GeomAPI_ProjectPointOnSurf` | `geometry::extrema::nearest_on_surface` |
| `STEPControl_Reader`, `STEPCAFControl_Reader` (names, colours, assembly structure) | `io::step::read`, `io::step::read_assembly` -> `io::step::Node` |
| `STEPControl_Writer`, `STEPCAFControl_Writer` | `io::step::write`, `io::step::write_assembly` |
| `BRepMesh_IncrementalMesh` | `meshing::mesh` |
| `BRepAlgoAPI_Fuse`, `Cut`, `Common` | `boolean::union`, `difference`, `intersection` (on solids, and on faces, shells, wires, edges and vertices with OCCT's rules for arguments of different dimensions) |
| `BRepAlgoAPI_Section` | `boolean::section` |
| `BOPAlgo_Builder` (the general fuse), `BOPAlgo_CellsBuilder`, `BOPAlgo_Splitter`; DRAW `baddobjects`, `baddtools`, `bfillds`, `bbuild`, `bsplit`, `bbop`, `bcbuild`, `bcadd`, `bcremove`, `bcremoveint` | `boolean::cells(model, shapes)` -> `Vec<Cell { inside, shape }>` (each piece with the shapes it lies inside; solid cells share the faces where they meet; a planar face across a solid splits it); YAML `partition`, `split` |
| `BOPAlgo_MakerVolume` (DRAW `mkvolume`, `-ni`) | `boolean::enclosed_solids(model, shapes)`, `enclosed_solids_unsplit` (one solid or a compound; internal shapes not kept) |
| `BRepOffsetAPI_ThruSections` (`AddWire`, `AddVertex`, `isSolid`, `ruled`); DRAW `thrusections` | `loft::loft(model, sections, &LoftOptions { solid, ruled, matching })` (sections: wires, or vertices at the ends; the first section repeated last closes the loft on itself; `Matching::ByRays` matches edges as `BRepFill_CompatibleWires` does, `AsGiven` as `CheckCompatibility(false)`); `Loft::generated` for `Generated`; YAML `loft` |
| `BRepOffsetAPI_MakeOffset` (planar wires, `GeomAbs_Arc`); DRAW `mkoffset` | `offset::planar_wire(model, wire, distance, Join::Arc)` (outwards for a positive distance from a closed wire, around an open one); YAML `offset_wire` (`distance`, or `length_fraction` of the wire's length) |
| `BRepOffsetAPI_MakeOffsetShape`, `BRepOffset_MakeOffset` (`GeomAbs_Intersection`; `GeomAbs_Arc`, rolling-ball faces where offsets part, except between faces moved by different distances); DRAW `offsetparameter`, `offsetload`, `offsetperform` | `offset::shape(model, shape, distance, join)` (the shape's topology kept: offsets that make faces collide or vanish are not handled); YAML `offset_shape` |
| `BRepOffsetAPI_MakeThickSolid` | `offset::thick_solid(model, solid, removed_faces, thickness, join)` (walls between the faces and their offsets, rims at the removed faces); YAML `thick_solid` |
| `BRepOffsetAPI_DraftAngle` (`Add` with a face, direction, angle, neutral plane), `Draft_Modification`; DRAW `depouille` | `draft::draft(model, shape, direction, &[DraftedFace { face, angle, neutral }])` (planes, and cylinders and cones about the pull direction; tangent neighbours that can be drafted follow; the topology is kept) |
| `BRepFeat_MakePrism`, `BRepFeat_MakeDPrism`, `BRepFeat_MakeRevol` (`Perform(length)`, `Perform(Until)`, `Perform(From, Until)`, `PerformUntilEnd`, `PerformThruAll`; fuse modes 0/1/2); DRAW `featprism`, `featdprism`, `featrevol`, `featperform`, `featperformval` | `feature::feature(model, base, profile, FeatureSweep::{Prism, DraftedPrism, Revolution}, FeatureMode::{Remove, Add, ToolOnly}, FeatureExtent::{Length, Until, FromUntil, UntilEnd, ThruAll})` (the tool trimmed by half-spaces at the limiting faces, then fused or cut; OCCT's local gluing and its selection of tool pieces are not reproduced) |
| `BRepFeat_MakeLinearForm`, `BRepFeat_MakeRevolutionForm` (ribs and slots); DRAW `featlf`, `featrf` | `feature::rib(model, base, wire, plane, RibSweep::{Linear, Revolution}, mode)` (the wire closed along the base's section on the side of OCCT's check point; a revolution turned as OCCT and cut to the slab) |
| `BRepOffsetAPI_NormalProjection`, `BRepAlgo_NormalProjection` (`Add`, `SetLimit`, `SetMaxDistance`, `SetParams` tolerance); DRAW `nproject` | `project::normal(model, shapes, target, &NormalProjection { tolerance, max_distance, within_faces })` (a compound of edges, one per piece on each face, with vertices of their own; projections collapsing to a point left out) |
| `Geom_OffsetSurface` (`Surface`, evaluation) | `geometry::surface::OffsetSurface` (`equivalent_surface`, `ParametricSurface`) |
| `BRepBuilderAPI_Sewing` (`SetTolerance`, `SetCutting`, `SetNonManifoldMode`, `Perform`, `SewedShape`); DRAW `sewing` | `sew::sew(model, shapes, &SewOptions { tolerance, cutting, non_manifold })` (manifold only; the result is a shell, a face, or a compound of those); YAML `sew` |
| DRAW `bfuseblend`, `bcutblend`; the general fuse builder (`BOPAlgo_Builder`) | YAML `filleted_union`, `filleted_difference`; `partition` |
| `BVH_Tree` and its builders (`BVH_BinnedBuilder`, ...) | `math::Bvh` (`new`, `with_options(BvhOptions { leaf_size, max_depth })`; `overlapping`, `containing`, `along_ray`, `overlapping_pairs`; `root`, `node_count`, `depth`, `leaf_sizes`, `sah_cost`) |
| `BVH_Traverse`, `BVH_PairTraverse` | `math::Selector` / `PairSelector` run by `Bvh::select` / `select_pairs` |
| `BVH_Distance`, `BVH_PairDistance` | `Bvh::nearest`, `Bvh::closest_pair` |
| `BVH_Box` (`Size`, `Center`, `Area`, `IsOut`, `Contains`, `Transformed`) | `math::BoundingBox` (`size`, `center`, `area`, `overlaps`, `contains_box`, `transformed` / `affine_transformed`) |
| `BVH_Tools` (`PointBoxSquareDistance`, `BoxBoxSquareDistance`, `PointBoxProjection`, `RayBoxIntersection`, `PointTriangleProjection`) | `BoundingBox::distance_to_point`, `distance_to_box`, `closest_point`, `ray_interval`; `math::closest_point_on_triangle` |
| `BRepAlgoAPI_BooleanOperation::SectionEdges` | `boolean::union_with_sections`, `difference_with_sections` (the result with its section edges) |
| `BRepFilletAPI_MakeFillet`, `MakeFillet::Add` (radius per edge) | `blend::fillet`, `blend::fillet_edges` |
| `BRepFilletAPI_MakeFillet::Add(UandR, E)` (evolving radius, `Law_Interpol`) | `blend::fillet_evolving` (`(u, r)` pairs per edge) |
| `BRepPrimAPI_MakeHalfSpace` | `primitives::half_space` (booleans treat an operand with an open shell as a half-space) |
| `BRepBuilderAPI_NurbsConvert` | `convert::to_nurbs` |
| `BRepBuilderAPI_Transform` (copying), `BRepBuilderAPI_GTransform`, `BRepTools_TrsfModification`, `BRepTools_GTrsfModification` | `modify::transform`, `modify::transform_affine` |
| `BRepBuilderAPI_MakeFace(wire)` | `make::planar_face` |
| `BRepBndLib::Add` + `Bnd_Box` | `properties::bounding_box` -> `math::BoundingBox` |
| `BRepBndLib::AddOptimal` | `properties::tight_bounding_box` |
| `ShapeUpgrade_UnifySameDomain` (`UnifyFaces`, `UnifyEdges`; DRAW `unifysamedom`, CadQuery `clean`) | `unify::unify_same_domain`, `unify_same_domain_with(UnifyOptions)` |
| `BRep_Tool::SameParameter`, `SameRange` | `Model::is_same_parameter`, `is_same_range` (per coedge) |
| `GProp_PrincipalProps::Moments` | `MassProperties::principal_moments` |
| `GC_MakeArcOfCircle` | `Circle::arc` |

## Operations

| OCCT | kentron |
|---|---|
| `Translated`, `Rotated`, `Scaled`, `Transformed` | `Transformable::translated`, `rotated`, `scaled`, `transformed` |
| `Mirrored(gp_Pnt)`, `Mirrored(gp_Ax1)`, `Mirrored(gp_Ax2)` | `mirrored_point`, `mirrored_axis`, `mirrored_plane(origin, normal)` |
| `Transform(...)`, `Mirror(...)` (in place) | `x = x.transformed(...)` |
| `IsEqual(other, tol)` | `approx_eq(other, tol)` |
| `IsNormal` | `is_perpendicular` |
| `Magnitude`, `SquareMagnitude`, `Modulus` | `length`, `length_squared` |
| `Crossed`, `CrossMagnitude`, `DotCross` | `cross`, `cross(..).length()`, `dot(..cross(..))` |
| `SetLinearForm(a, u, b, v)` | `a * u + b * v` |
| `Reversed()` | unary `-` |
| `gp_Vec(p1, p2)` | `p2 - p1` |
| `BaryCenter` | `Point::barycenter(&[(p, w), ...])` |
| `AngleWithRef` | `Direction::signed_angle` |
| `SetDisplacement(from, to)` | `Transform::frame_to_frame(&from, &to)` |
| `Coefficients(a, b, c, d)` out-params | tuple or array return (`Plane::equation`, `Torus::implicit_coefficients`) |
| Constructor throwing `Standard_ConstructionError` | `try_*` / `Option`-returning constructor |
