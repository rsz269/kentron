# CadQuery scripts found failing, diverging or slow on pykentron

Scripts from public CadQuery projects on GitHub (2026-10-18 and 2026-10-19) that pykentron
fails on, builds differently from CadQuery 2.8, or builds substantially slower. Each
project's scripts keep its licence, copied beside them (`<project>/LICENSE`). Scripts under
other licences than MIT, Apache, BSD and the Unlicense, or none, are not copied: their
cases are only listed.

The source of every case and what went wrong is in `cases.tsv`:

- the project, its repository and commit;
- the licence;
- the kind: `fails`, `diverges` or `slower`;
- the suspected area;
- both kernels' wall times;
- the first differences or the error.

## How the cases were found

`harness/` holds the tools.

- `runner.py` runs a script as `__main__` in a fresh, isolated process, on pykentron's shim or on CadQuery. It records what the script shows (`show_object`, `show`), exports, or leaves at module level: volume, area, bounding box, validity, and solid and face counts. The volume is its solids' (OCCT gives wires and faces in a compound volumes of their own); CadQuery's bounding box is OCCT's tight one (its default is grown by tolerances and control polygons); objects shown unnamed are numbered among themselves, exports not counted (CadQuery records some exports twice). Since 2026-10-23.
- `hunt.py` runs up to 30 scripts per project on both kernels, with a 120 s limit, and compares the results.
- `analyze.py` takes the first 20 projects in which CadQuery itself built anything, and gives each flagged case a suspected area.

A case is flagged if any of these holds:

- pykentron raises where CadQuery doesn't;
- a volume or area differs by more than 0.1%, or a bounding box by more than 1% of its size;
- the solid count differs;
- the result is invalid where CadQuery's is valid;
- pykentron is more than 2× and 0.5 s slower.

Scripts using unseeded randomness are left out, since no two runs of them agree. So are scripts that call OCCT's Python bindings (OCP) for more than the shapes CadQuery hands them: cadquery-massembly's examples (`OCP.Geom` curves, 5 scripts) and Cadquery_Experiments' BarbFitting (`OCP.gp.gp_Quaternion` between two vectors), dropped 2026-10-23. They test OCCT's API, not CadQuery's. LatticeQuery's varying_Ls draws its sizes from `random()` and is left out too.

Many scripts import their project's own library: cq_warehouse, cqterrain, cadqueryhelper, cqindustry, cq_gears, cqMore and cq-electronics, and from the second round LatticeQuery and gridfinity-cadquery's package. To run those, install the library without its dependencies (`--no-deps`), so that it can't replace CadQuery.

A second round (2026-10-19) took 57 more projects from GitHub's repository search (`cadquery`, `topic:cadquery`, and `cadquery` with gridfinity, keyboard, enclosure, printable, parametric and the like), leaving out lists, AI pipelines and datasets: 382 scripts, of which CadQuery built 187. pykentron built 129 of those as CadQuery does. Of the other 58, the 53 cases at the end of `cases.tsv`; left out are two scripts writing to a directory that does not exist, one reading a file on its author's disk, one using unseeded randomness and one showing nothing to compare.

A third round (2026-10-20) took 10 more projects from the search (`cadquery` with 3d print, parametric, stl, design and the like): 61 scripts, of which CadQuery built 44. pykentron built 29 of those as CadQuery does; 14 of the others are the last cases in `cases.tsv` (one more imports `random`, unused, and is left out). For these, the runner stands in for the VS Code viewer's module (`ocp_vscode`) as it does for CQ-editor's `show_object`.

## Scripts that build as on CadQuery

`passing.tsv` lists every script of all three rounds that CadQuery built and pykentron builds as it does (the same checks as for the cases, run again on 2026-10-20): 274 scripts, 189 of them copied here under their projects' licences, with the project's own modules they import. Each row gives both kernels' wall and CPU times (CadQuery's from the hunt, in a fresh process each; it runs several threads, so its CPU time is often over its wall time). In all, pykentron took 194 s where CadQuery took 469 s (CPU: 284 s against 1408 s), and was faster on 238 of the 274.

`harness/timing.py` runs the copied ones on pykentron again, each in a fresh process, and prints their times against CadQuery's, with the totals: `python3 harness/timing.py --python PYTHON [FILTER...] [--repeat N]`. `PYTHON` needs the libraries some scripts import (cq_warehouse, cqterrain and the like, installed as above, and numpy); CadQuery itself is not used.

## Since

Fixed since the cases were found (running them again):

- OCP names cq_warehouse calls on shapes and vectors (`BRepBuilderAPI_Copy`,
  `ShapeUpgrade_UnifySameDomain`, `TopoDS_Shape`'s and `gp_Vec`'s methods, ...): its
  sprocket and chain examples now build as on CadQuery.
- Smoothed spline approximation (what `parametricCurve` makes): cq_warehouse's threads and
  the contrib cylindrical gear now build as on CadQuery.
- A wire's point at a fraction of its length, at its end: was searched for past it (2 s a
  point on a helix); cq_warehouse's fastener examples took 34–102 s, now about 1 s.
- Sewing faces whose edges collapse to a point there (a ruled face between curves that
  meet): the faces lost their area. The contrib thread's solids are now right (CadQuery's
  own volume of them is 3% off: 13.35 where the surfaces enclose 12.99).
- Fillets and chamfers along an arc whose circle comes back into the solid: the tool ran
  round into it.
- Fillets and chamfers given edges where the faces meet smoothly (a hole's seam, picked
  by `edges("|Z")`): those left out, as OCCT leaves them (CadQuery's Ex003, the Raspberry
  Pi board).
- Fillets and chamfers on several solids, one shape placed twice among them: each solid on
  its own (cqterrain's rooms).
- A surface of revolution whose arc lies across its axis from the arc's centre (a concave
  fillet round a thin shaft): points on it were found off it, and booleans cut it wrongly
  (the contrib resin mould, the svg-export mould).
- Wires joined from mirrored halves: the halves made to wind one way, the edge out and back
  along the mirror line dropped (cqterrain's pointed arches); a wire given twice to a face
  counted once (cq_warehouse's domed nuts); free wires transformed (mirrored).

- Lofts matched as OCCT matches them (by rays from the sections' centres where their edge
  counts differ); closed sections in planes near square to each other turned and started
  by their vertices' distances, not by their normals or the nearest start (CadQuery's
  Ex018, the chute loft).
- Sweeps along spines with corners mitred or rounded as CadQuery's `transition` says (it
  was ignored: a closed spine's sweep came out empty), and carried in the spine's Frenet
  frame where `isFrenet` (CadQuery's Ex023 and Ex025, both Tetrakaidecahedron scripts, now
  15–20 times faster than on CadQuery).
- A planar face's centre exactly in its plane: faces in one plane, sorted by height, keep
  their order (CadQuery's Ex026 selected the lid's face for the rim's).

- Shapes translated or rotated, then `located`: the location placing the moved shape, as
  OCCT's moves the copy's geometry first (cadqueryhelper's rotated and scheme grids,
  cqterrain's books and bookcases).
- Tapered extrusions of outlines whose straight sides are split in line (cadqueryhelper's
  cell grids).

- Spline surfaces approximating grids of points (`Face.makeSplineApprox`), and OCP's sewing
  (`BRepBuilderAPI_Sewing`) as cq_gears calls it: its planetary gearsets build as on
  CadQuery; single spur, helical and bored gears too.

- Sweeps through several sections: blended smoothly through them (a parabola through three,
  as OCCT's), not linearly between neighbours; sections of different numbers of edges
  matched as a loft's; the end caps facing out (CadQuery's Ex024 but its last sweep, whose
  sections lie off the one edge it passes as the path; ten more of OCCT's cases pass).
- Shape constructors called on subclasses make the subclass; vectors without a length (numpy
  keeps arrays of them as objects); planes of negative size: cqMore's dice now build (their
  engraved text differs, a font's), and the spindle vacuum shroud (CadQuery's result of it
  is empty).

- Fusing many shapes at once: pairwise, then pairs of pairs, not each added in turn to a
  growing whole (cqindustry's bridge, 263 pieces: 76 s, now 3 s, as on CadQuery).

- Fillets and chamfers of straight edges meeting at a corner concave on the face they share
  (a pocket's or opening's rim): mitred there, as OCCT's, not left square with small faces
  filling the corner (the contrib shelled cube).

- A hole through a solid's side whose edge it crosses twice, once exactly at a point the
  edge is sampled at: the other crossing was missed, and the hole not cut (the contrib 3D
  printer extruder support; what still differs there is a workplane on three coplanar faces
  facing two ways, taking the first's normal: OCCT's face order, not kentron's).
- Fillets of the mitre where two fillets meet: carried on from the edge going on smoothly
  into it, narrowing to nothing where its faces turn tangent; balls nearly as large as the
  mitre's curvature taken in their own planes across (cadquery-models' lid fixer: CadQuery's
  own fillet there strays 0.31 from the ball, 0.23% of the volume).
- Fillets running from a straight side round a rounded corner, joined across the one ball
  touching both (they were joined loosely, their edges 0.5 apart); faces whose boundary
  edges meet only within their tolerances split by booleans (cadquery-models' foot block).
- Intersections with several shapes, with all of them (OCCT's common), not with each in
  turn (cqindustry's round platform and can tower stairs).

- Projections, fillings and thickening, as cq_warehouse calls them through OCP: wires
  projected onto shapes (along a direction or from a point), faces spanning loops of edges
  (the plate of least bending energy through them and given points), holes cut in curved
  faces, curved faces thickened, lines crossing shapes' faces (cq_warehouse's cylinder with
  holes now builds as on CadQuery; its heatset nuts too). Swept faces of curves with no
  tangent at an end (font outlines') no longer meet other surfaces all along a sphere's seam.
  CadQuery's own fillings stray from their points (by up to 1.4 in its interpolated plates),
  which kentron's pass through: Ex101's plates differ by 0.2–0.6%; fillings approximating a
  surface from its outline (cq_warehouse's projected flag faces) by about 1%; a filling's
  normal, and so the way it thickens, follows the loop in kentron (cq_warehouse's squares
  mapped on a sphere thicken the other way). OCCT's follows no rule to match: it ignores
  which way the loop runs, and as a square loop tilts about an axis in its plane its normal
  flips back and forth (every 15° in places), at some tilts with the side of the loop's
  plane its interior point is on. Its projection example's conical
  case projects text from a point in the text's own plane (every projected outline folded
  onto one circle): it stops there.

- Spline surfaces fitted to grids with few points one way (no more than a span of the
  greatest degree has poles): through them exactly that way, as OCCT's, so that gear teeth's
  faces fitted to grids sharing a row meet along it; a shape split by something not cutting
  it is itself, not a compound (as OCCT's); OCP's surface intersections (`GeomAPI_IntSS`, of
  faces), edges and wires made from them: cq_gears' worm gear and gear rack build as on
  CadQuery.

- Lofts whose end sections lie in no one plane (wires projected onto a curved face) capped
  by fillings on their edges: cq_warehouse's teapot handle lofts (the teapot then stops at
  the spout's union with the pot, faces on one surface); fillets whose sections lose their
  arc, and sections of surfaces all but meeting at a point, fail with an error, not a crash
  (cq_warehouse's bottle funnel: a fillet as large as its outline's corners and curves, still
  open).

- Speed: a fillet whose ball only grazes the faces beyond its own rolled on its own faces
  first (rolling it on over them, the slower way, only where that fails): cadquery-models'
  foot block 3.4 s, now 0.5 s; an ellipse of equal radii made a circle, so that what is
  extruded or offset from it stays a cylinder or circle (cqindustry's round platform 5.6 s,
  now 0.25 s; its can tower stairs 5 s, now 2 s); sweeps through several sections as compact
  along the spine as the sections' blend allows (spindle's chimney 1.2 s, now 0.9 s).

Three scripts build differently on CadQuery for reasons of its own: cqterrain's magnet
template queues its pending wire again (`toPending()`), and CadQuery's face counts it twice
(double the volume); cqterrain's buildings stack storeys as touching solids, which CadQuery
leaves 14 to 22 solids and kentron's union one; the contrib drawer frame fuses dovetails
into its walls with `glue=True`, which assumes shapes only touch, and they overlap.

With these, five more of OCCT's cases pass (four booleans and a Tcl bug script).

Still open among them: the thread solids' unions with their core (nearly coincident
surfaces, 0.001 apart) and two meshing helical gears' union (teeth touching), cq_warehouse's bolts (their threads' unions too, by the look of
them), the contrib remote enclosure (a countersink through a wall thinner than it is deep:
the cut loses the countersink's cone), and the rest of `cases.tsv`. 45 of the 85 now build
as on CadQuery.
