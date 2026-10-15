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
drafts, offsets and thick solids, sewing and healing, text from TrueType and OpenType
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

## Tested against

- **OCCT's DRAW test scripts:** 13,074 of its 13,183 scripts (the rest exercise viewers and
  data exchange). 9,596 are converted by `tools/draw2rs.py` into Rust tests and YAML cases;
  the 3,478 that use Tcl (conditions, loops, procedures) run as they are, on a Tcl
  interpreter. Many need OCCT's external test dataset (not in this repo).
- **OCCT's GTests:** the computational-geometry and topology toolkits, translated by hand.
- **BRL-CAD:** NIST's MBE PMI test parts and other STEP files from its source tree, with
  OCCT's reading of them as the reference, plus its `libbg` geometry tests.
- **CadQuery's own test suite** (v2.8.0), run unchanged against pykentron: 547 of 727 tests
  pass. Most of the rest exercise OCCT's document formats (OCAF/XDE), exporters kentron
  lacks (glTF, VRML, VTK), or OCCT-specific behaviour.
- **44 CadQuery models from a private project** build on pykentron; for the 42 that return
  their parts, every part is valid and within 0.1% of CadQuery's volume.

About 10,300 tests run (with OCCT's dataset). About 1,700 of them are known gaps, marked as
expected failures (`crates/occt-tests/expected_failures.txt` for OCCT's scripts): they pass
while they fail, so `cargo test` is green, and fail once they pass, so a fix is noticed
as surely as a regression.

## Performance

A complex five-part CadQuery model, built on pykentron and on CadQuery 2.8 (OCCT):

| | wall time | CPU time | peak memory |
|---|---|---|---|
| pykentron (`fast` profile) | 4.9 s | 5.9 s | 460 MB |
| CadQuery 2.8 / OCCT | 41.8 s | 122.4 s | 1,097 MB |
| | **8.5× faster** | **21× less** | **2.4× less** |

## Memory or speed: profiles

pykentron frees shapes that Python no longer refers to. How much it frees, and which
allocator it uses, is chosen with `PYKENTRON_PROFILE`:

| profile | what it does | the same model, against no collection |
|---|---|---|
| `low-memory` | also drops cached data of removed shapes | 58% of the peak memory, ~3% slower |
| `balanced` (default) | removes unreachable shapes | 68% of the peak memory, ~2% slower |
| `fast` | as `balanced`, with mimalloc as allocator | 89% of the peak memory, ~14% faster |

Rust programs using kentron directly make the same choice by declaring their allocator and
passing `CollectOptions` to `Model::collect_garbage_with`.

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
