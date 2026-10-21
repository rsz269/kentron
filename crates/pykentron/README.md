# pykentron

Python bindings for kentron, with a CadQuery-compatible API, so that CadQuery scripts run
on kentron unchanged.

## Installing

```
pip install pykentron              # the kernel and pykentron's own CadQuery-like API
pip install "pykentron[cadquery]"  # also what CadQuery's own layers need, to run CadQuery scripts
```

Python 3.10 or later. The wheels on PyPI are built for x86-64-v2 Linux machines (Intel from 2008, AMD from 2011; build from source with an empty `RUSTFLAGS=` for older ones).

## Using it

```python
import pykentron as cq
box = cq.Workplane("XY").box(10, 20, 30).faces(">Z").workplane().hole(5)
print(box.val().Volume())
```

To run CadQuery scripts unchanged, import `pykentron.shim` first: `import cadquery` then gives
CadQuery's own pure-Python layers (`Workplane`, `Sketch`, selectors, `cqgi`, free functions,
assemblies: copied unchanged from CadQuery v2.8.0, Apache-2.0, see
`pykentron/cadquery_shim/cadquery/VENDORED.txt`) on pykentron's shapes. Or run a script as it
is with `pykentron-run script.py`. Nothing named `cadquery` is installed, so a real CadQuery
beside pykentron is unaffected unless the shim is activated (the two cannot share one
process).

Code written against OCP, as libraries like cq_warehouse are, finds some of OCCT's names:

- `Shape.wrapped` answers `TopoDS_Shape`'s `Moved`, `Located`, `Location` (also placing the
  shape in place), `Orientation`, `Oriented`, `Reversed`, `ShapeType`, `IsSame`, `IsEqual`,
  `IsPartner`, `IsNull`;
- `Vector.wrapped` answers `gp_Vec`'s methods, `Location.wrapped` `IsEqual` and
  `Transformation`, and `Matrix.wrapped` `Trsf()`;
- the classes `BRepBuilderAPI_Copy`, `BRepBuilderAPI_Transform`, `BRepBuilderAPI_Sewing`,
  `BRepBuilderAPI_MakeFace` (a face from a planar wire, or a face with holes added),
  `ShapeUpgrade_UnifySameDomain`, `Bnd_Box` with `BRepBndLib.AddOptimal_s`,
  `BRepTools.UVBounds_s` and `Clean_s`, the empty `TopoDS_Shape()`, the `TopAbs` numbers and
  the `gp` classes;
- projections and fillings: `BRepProj_Projection` (wires onto shapes, along a direction or
  from a point), `BRepOffsetAPI_MakeFilling` (faces spanning a loop of edges through points:
  the plate of least bending energy, as `Face.makeNSidedSurface` and `interpPlate` make),
  `BRepOffset_MakeOffset` (thickening a face into a solid), `BRepIntCurveSurface_Inter` with
  `gce_MakeLin` (where a line crosses a shape's faces), `ShapeAnalysis_FreeBounds`'
  `ConnectEdgesToWires_s` with `TopTools_HSequenceOfShape`, `GeomAPI_IntSS` (on the faces
  `BRepAdaptor_Surface` stands for), `BRepBuilderAPI_MakeEdge` and `MakeWire`, `TopoDS`'
  casts, `ShapeFix_Face` (keeping the face), and `StdFail_NotDone` (as `ValueError`, what
  kentron's failures raise);
- `occ_impl.shapes`' `downcast`, `fix` and `shapetype`.

Other names, of pykentron or OCP, import but raise `NotImplementedError` when used.

## Building from source

```
pip install ./crates/pykentron                 # with maturin, as the PyPI wheels are built
python3 crates/pykentron/build.py              # in place, optimized (LTO), into python/pykentron/_kentron.abi3.so
python3 crates/pykentron/build.py --quick      # in place, the tests' release profile
```

The in-place build installs the extension by renaming, so interpreters that have it loaded
keep running; put `crates/pykentron/python` on `sys.path` to use it. The CadQuery layers need
`pyparsing`, `multimethod`, `runtype`, `typing_extensions` and `numpy` (and `nlopt`, `casadi`
for constraints), as CadQuery does: the `cadquery` extra.

## Memory and speed

Shapes no Python object refers to any more are collected from time to time. How much is
freed, and which allocator kentron's memory comes from, is a profile, chosen by
`PYKENTRON_PROFILE` before `pykentron` is imported:

| profile | allocator | collection | a complex five-part CadQuery model, against `balanced` |
|---|---|---|---|
| `low-memory` | the system's | also forgets cached values of removed shapes | 88% of the peak memory, about as fast |
| `balanced` | the system's | removes unreachable shapes | (360 MB at the peak, 3.9 s) |
| `fast` (default) | mimalloc | as `balanced` | 132% of the peak memory, ~11% faster |

`pykentron.profile()` gives the profile in use and `pykentron.set_profile(name)` switches
between `low-memory` and `balanced` at any time; the allocator is chosen once, at startup,
so `fast` is only had through the variable. Building without the `mimalloc` feature
(`cargo build -p pykentron --no-default-features`) drops mimalloc and the `fast` profile
(`balanced` is then the default).

## What is there

- `Workplane`: sketching (polylines, lines, arcs, rectangles, circles, polygons, points,
  workplanes and transforms), `extrude`/`cutBlind`/`revolve`/`box`/`cylinder`,
  booleans, `fillet`, selection (`faces`, `edges`, ... with string selectors `|`, `#`,
  `+`, `-`, `<`, `>`, `%TYPE`), placement.
- Shapes (`Vertex`, `Edge`, `Wire`, `Face`, `Shell`, `Solid`, `Compound`): construction
  (boxes, cylinders, cones, spheres, tori, prisms, revolutions, lofts, `offset2D`,
  shell offsets and thick solids), booleans, `clean`, transforms and mirrors, queries
  (volume, area, length, centres, bounding boxes, validity, normals, geometry types,
  `distance`, `isInside`, `tessellate`). As in CadQuery, a shape's bounding box is tight
  until it is tessellated, then bounded by its meshes, widened by their deviation.
- Files: STEP import and export (shapes, or an `Assembly` with names, colours and
  placements), BRep import and export, STL export.
- Splines: interpolated, approximated (`makeSplineApprox`, smoothed too, as
  `parametricCurve` makes them), Bézier.
- `Assembly`, `Location`, `Plane`, `Vector`, `Color` (by OCCT's colour names, which are
  X11's), `BoundBox`.

Not yet: DXF import, SVG and other exporters, `cadquery.vis.show` (warns).

## Checking it

- `tests/test_against_cadquery.py` compares pykentron with the real CadQuery case by case
  (run it with a Python that has CadQuery: `PYTHON
  crates/pykentron/tests/test_against_cadquery.py`).
- `scripts/cadquery_build.py --project DIR MODULE` runs a CadQuery project's `build()` on
  pykentron and reports each part's volume, faces and validity; `--oracle` runs it on
  CadQuery instead.
- `scripts/cadquery_sweep.py --project DIR --python PYTHON` builds every module of a
  project both ways and checks that each part is valid and within 0.1% of CadQuery's
  volume (`--oracle-cache FILE` keeps CadQuery's results, which do not change).
- `scripts/cadquery_tests.py` runs CadQuery's own test suite (release v2.8.0, fetched from
  GitHub into `target/cadquery-tests` with a virtual environment of its own, made with `uv`)
  against pykentron, unchanged: `python3 crates/pykentron/scripts/cadquery_tests.py [pytest
  arguments]`. It fails on any test failing that is not in `tests/cadquery_baseline.txt` and
  reports those of the baseline that pass; `--update-baseline` records the current failures.
  The suite reaches CadQuery's internals (`cadquery.occ_impl.*`, `cadquery.func`, OCP): the
  shim resolves those names to pykentron's where it has them, and makes the rest placeholders
  that raise NotImplementedError when used.
- `tests/found/` holds CadQuery scripts from public projects that pykentron failed on,
  built differently from CadQuery, or built much more slowly, with where each came from and
  what went wrong (`cases.tsv`), and the tools that found them (`harness/`: each script run
  on both, in fresh isolated processes, and compared).
