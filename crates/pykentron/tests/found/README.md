# CadQuery scripts found failing, diverging or slow on pykentron

Scripts from public CadQuery projects on GitHub (2026-10-18) that pykentron fails on,
builds differently from CadQuery 2.8, or builds substantially slower. Each project's scripts
keep its licence, copied beside them (`<project>/LICENSE`).

The source of every case and what went wrong is in `cases.tsv`:

- the project, its repository and commit;
- the licence;
- the kind: `fails`, `diverges` or `slower`;
- the suspected area;
- both kernels' wall times;
- the first differences or the error.

## How the cases were found

`harness/` holds the tools.

- `runner.py` runs a script as `__main__` in a fresh, isolated process, on pykentron's shim or on CadQuery. It records what the script shows (`show_object`, `show`), exports, or leaves at module level: volume, area, bounding box, validity, and solid and face counts.
- `hunt.py` runs up to 30 scripts per project on both kernels, with a 120 s limit, and compares the results.
- `analyze.py` takes the first 20 projects in which CadQuery itself built anything, and gives each flagged case a suspected area.

A case is flagged if any of these holds:

- pykentron raises where CadQuery doesn't;
- a volume or area differs by more than 0.1%, or a bounding box by more than 1% of its size;
- the solid count differs;
- the result is invalid where CadQuery's is valid;
- pykentron is more than 2× and 0.5 s slower.

Scripts using unseeded randomness are left out, since no two runs of them agree.

Many scripts import their project's own library: cq_warehouse, cqterrain, cadqueryhelper, cqindustry, cq_gears, cqMore and cq-electronics. To run those, install the library without its dependencies (`--no-deps`), so that it can't replace CadQuery.

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

With these, five more of OCCT's cases pass (four booleans and a Tcl bug script).

Still open among them: the thread solids' unions with their core (nearly coincident
surfaces, 0.001 apart), cq_warehouse's bolts (their threads' unions too, by the look of
them), the contrib remote enclosure (a countersink through a wall thinner than it is deep:
the cut loses the countersink's cone), and the rest of `cases.tsv`. 20 of the 85 now build
as on CadQuery.
