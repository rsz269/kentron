# kentron

[![CI](https://github.com/rsz269/kentron/actions/workflows/ci.yml/badge.svg)](https://github.com/rsz269/kentron/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/kentron.svg)](https://crates.io/crates/kentron)
[![PyPI](https://img.shields.io/pypi/v/pykentron.svg)](https://pypi.org/project/pykentron/)

A boundary-representation CAD kernel written from scratch in Rust, with Python bindings that
aim to run CadQuery scripts unchanged.

## Why

Because I got tired of my complex-ish CadQuery scripts taking several minutes to finish and
eating up 1 GB+ of memory, so I figured why not spin up my own kernel.

## What it is

kentron models solids the way Open CASCADE Technology (OCCT) does: as faces bounded by
edges and vertices, over exact curves and surfaces (analytic and NURBS). It covers
primitives, extrusions and revolutions, lofts and pipes, booleans, fillets, chamfers,
drafts, offsets and thick solids, removing features, sewing and healing, meshing, text from
TrueType and OpenType
fonts, and BRep, STEP and STL files. All of
it is kentron's own code: OCCT is the reference to match, not a dependency.

- **`crates/kentron`**: the kernel, [on crates.io](https://crates.io/crates/kentron). Pure
  Rust, with two optional dependencies: `rayon`, for running independent steps in parallel,
  and `ttf-parser`, for text (the `text` feature).
- **`crates/pykentron`**: Python bindings with a CadQuery-compatible API,
  [on PyPI](https://pypi.org/project/pykentron/). CadQuery scripts run on kentron by
  importing `pykentron` in place of `cadquery`, or unchanged with `pykentron-run`.
- **`crates/occt-tests`**: OCCT's own tests, converted to run against kentron.

## Aims and scope

- **OCCT's functionality, not its API.** Names, conventions and data structures are
  kentron's own. All topology lives in one arena (`topology::Model`) addressed by typed ids,
  with no reference counting. `docs/occt-mapping.md` maps OCCT concepts to kentron's.
- **Matching OCCT where it matters:** the same shapes, volumes and validity on OCCT's test
  cases. Where OCCT's expected values are themselves wrong, the converted tests say so and
  why.
- **Running real CadQuery code** faster and in less memory than CadQuery on OCCT.
- **Currently out of scope:** OCCT's visualization, hidden-line removal, OCAF/XDE documents
  and its internal bookkeeping.

It is pre-1.0. Some OCCT features are still missing, for example sweeps with several
sections of different shapes, offsets that change a shape's topology, evolved solids, and
fillets meeting in some kinds of corner (four or more faces at a vertex, convex and concave
fillets together). Some rare pathological
cases may still run very slowly or fail to converge.

## How it was made

kentron is 99% vibe coded: written almost entirely by an AI coding agent (Claude, in Claude
Code), with a human setting the direction and judging the results. The tests are what keep
it honest. Every change is checked against the test sources below before it is committed,
and no change may make a passing test fail.

<!-- effort -->
So far (2026-09-25 to 2026-10-18): about 104 hours of the agent's active work, 8.5 million
tokens written by it, and 4.1 billion tokens processed in all (nearly all of them its
conversation's context, read back from the cache).
<!-- /effort -->

## Tested against

- **OCCT's DRAW test scripts:** 16,311 of its 16,420 scripts in the modelling groups (the
  rest exercise viewers and data exchange). 10,245 are converted by `tools/draw2rs.py` into
  Rust tests and YAML cases; the 6,066 that use Tcl (conditions, loops, procedures) run as
  they are, on a Tcl interpreter. Many need OCCT's external test dataset (not in this
  repo).
- **OCCT's GTests:** the computational-geometry and topology toolkits, translated by hand.
- **BRL-CAD:** NIST's MBE PMI test parts and other STEP files from its source tree, with
  OCCT's reading of them as the reference, plus its `libbg` geometry tests.
- **CadQuery's own test suite** (v2.8.0), run unchanged against pykentron: 558 of 727 tests
  pass. Most of the rest exercise OCCT's document formats (OCAF/XDE), exporters kentron
  lacks (glTF, VRML, VTK), or OCCT-specific behaviour.
- **44 CadQuery models from a private project** build on pykentron; for the 42 that return
  their parts, every part is valid and within 0.1% of CadQuery's volume.
- **CadQuery scripts from 20 public projects** (up to 30 each), run on both: the 85 that
  pykentron failed on, built differently or built much more slowly are kept, with where each
  came from and what went wrong, in `crates/pykentron/tests/found/`, the cases still to fix.

About 12,700 tests run (with OCCT's dataset). About 1,500 of them are known gaps, marked as
expected failures (`crates/occt-tests/expected_failures.txt` for OCCT's scripts): they pass
while they fail, so `cargo test` is green, and fail once they pass, so a fix is noticed
as surely as a regression.

## Performance

A complex five-part CadQuery model, built on pykentron under each of its profiles and on
CadQuery 2.8 (OCCT), median of 5 runs:

| | wall time | CPU time | peak memory |
|---|---|---|---|
| pykentron, `fast` profile | 4.8 s | 5.9 s | 420 MB |
| pykentron, `balanced` profile (default) | 5.4 s | 6.4 s | 339 MB |
| pykentron, `low-memory` profile | 5.4 s | 6.4 s | 292 MB |
| CadQuery 2.8 / OCCT | 40.9 s | 125.8 s | 1,094 MB |

With the `fast` profile, that is 8.6× faster than CadQuery, with 21× less CPU time and 2.6×
less memory.

pykentron frees shapes that Python no longer refers to. How much it frees, and which
allocator it uses, is chosen with `PYKENTRON_PROFILE`:

- `low-memory`: also drops cached data of removed shapes.
- `balanced` (default): removes unreachable shapes.
- `fast`: as `balanced`, with mimalloc as allocator.

Rust programs using kentron directly make the same choice by declaring their allocator and
passing `CollectOptions` to `Model::collect_garbage_with`.

## Tracing sections of free-form surfaces

Where a boolean operation meets a B-spline, Bezier or offset surface, the curves where the
faces meet are traced numerically. `KENTRON_MARCH` chooses how; it is read once per process:

| value | how sections are traced |
|---|---|
| unset, `1`, `on` (default) | by marching along both surfaces at once from seeds found on a coarse grid, falling back to the grid where a branch cannot be followed |
| `0`, `off` | on a grid over one surface, where the other's signed distance changes sign |

Marching is several times faster on large free-form faces (a boolean from the OCCT suite,
`bugs/modalg_7/bug29038`: 16.6 s on the grid, 4.0 s marching); the grid is kept for
comparison.

## Installing

- Rust: `cargo add kentron`.
- Python 3.10 or later: `pip install pykentron`, or `pip install "pykentron[cadquery]"` to
  run CadQuery scripts (it adds what CadQuery's own Python layers need). Wheels are built
  for Linux on x86-64; elsewhere pip builds from the source distribution, which needs Rust.

```python
import pykentron as cq
box = cq.Workplane().box(10, 20, 30).edges("|Z").fillet(2)
```

An existing CadQuery script runs unchanged with `pykentron-run script.py`, or after
`import pykentron.shim` (before anything imports `cadquery`).

## Building and testing

- Rust 1.88 or later.
- `cargo test` runs everything. Most of OCCT's converted tests read its test dataset, an
  external download: see `crates/occt-tests/README.md`. Without it those tests fail.
- `tools/run_tests.py` runs the same tests with every test binary at once (about 20 s on 16
  threads, against about 30 s for `cargo test`, which runs them one after another).
  `--quick` leaves out the tests expected to fail (about 12 s); a fix that makes one pass
  then goes unnoticed.
- pykentron from source: `python3 crates/pykentron/build.py` (see
  `crates/pykentron/README.md`).

## More

- `crates/pykentron/README.md`: building and using the Python bindings.
- `crates/occt-tests/README.md`: running OCCT's tests, and getting their data.

## Licence

kentron and pykentron are dual-licensed under the [MIT](LICENSE-MIT) and
[Apache 2.0](LICENSE-APACHE) licences, at your option. Two parts keep the licences of the
code they come from:
- `crates/occt-tests`, translated from OCCT's and BRL-CAD's tests: LGPL-2.1 (OCCT's with the
  Open CASCADE exception). Nothing else depends on it.
- The CadQuery code copied into pykentron
  (`crates/pykentron/python/pykentron/cadquery_shim/cadquery` and
  `crates/pykentron/python/pykentron/occ_impl/exporters`, each with its
  `LICENSE.CadQuery` and a list of the files copied in `VENDORED.txt`): Apache-2.0.
