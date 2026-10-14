# pykentron

Python bindings for kentron, with a CadQuery-compatible API, so that CadQuery scripts run
on kentron unchanged.

## Installing

```
pip install pykentron              # the kernel and pykentron's own CadQuery-like API
pip install "pykentron[cadquery]"  # also what CadQuery's own layers need, to run CadQuery scripts
```

Python 3.10 or later. The wheels on PyPI are built for any x86-64 Linux machine.

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

Names pykentron does not have yet, and the OCP names scripts use directly (beyond `Bnd_Box`,
`BRepBndLib.AddOptimal_s`, `BRepTools.Clean_s`), import but raise `NotImplementedError` when
used.

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
