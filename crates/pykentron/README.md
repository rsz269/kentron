# pykentron

Python bindings for kentron, with a CadQuery-compatible API, so that CadQuery scripts run
on kentron unchanged.

## Building

```
python3 crates/pykentron/build.py             # optimized (LTO, AVX2/FMA) into python/pykentron/_kentron.abi3.so
python3 crates/pykentron/build.py --portable  # for any x86-64 CPU
python3 crates/pykentron/build.py --debug     # unoptimized
```

The extension is abi3 (Python 3.9 and later); the build installs it by renaming, so
interpreters that have it loaded keep running.

## Using it

Put `crates/pykentron/python` on `sys.path` and `import pykentron as cq`. To run a CadQuery
script without editing it, also put `crates/pykentron/python/cadquery_shim` first on the
path. Its `cadquery` package is CadQuery's own pure-Python layers (`Workplane`, `Sketch`,
selectors, `cqgi`, free functions, assemblies: copied unchanged from CadQuery v2.8.0 by
`scripts/vendor_cadquery.py`, Apache-2.0, see `cadquery_shim/cadquery/VENDORED.txt`) on
CadQuery's `occ_impl` layer taken from pykentron (`pykentron.occ_impl`, `pykentron.shapes`,
`pykentron.geom`). Names pykentron does not have yet, and the `OCP` names scripts and tests
use directly (beyond `Bnd_Box`, `BRepBndLib.AddOptimal_s`, `BRepTools.Clean_s`), import but
raise `NotImplementedError` when used. The layers need `pyparsing`, `multimethod`,
`runtype`, `typing_extensions` and `numpy` (and `nlopt`, `casadi` for constraints), as
CadQuery does.

## Memory and speed

Shapes no Python object refers to any more are collected from time to time. How much is
freed, and which allocator kentron's memory comes from, is a profile, chosen by
`PYKENTRON_PROFILE` before `pykentron` is imported:

| profile | allocator | collection | a complex five-part CadQuery model, against no collection |
|---|---|---|---|
| `low-memory` | the system's | also forgets cached values of removed shapes | 62% of the peak memory, ~5–10% slower |
| `balanced` (default) | the system's | removes unreachable shapes | 72% of the peak memory, ~3% slower |
| `fast` | mimalloc | as `balanced` | the same peak memory, ~17% faster |

`pykentron.profile()` gives the profile in use and `pykentron.set_profile(name)` switches
between `low-memory` and `balanced` at any time; the allocator is chosen once, at startup,
so `fast` is only had through the variable. Building without the `mimalloc` feature
(`cargo build -p pykentron --no-default-features`) drops mimalloc and the `fast` profile.

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
- `Assembly`, `Location`, `Plane`, `Vector`, `Color`, `BoundBox`.

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
