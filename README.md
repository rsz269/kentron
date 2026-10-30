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

- **[`crates/kentron`](crates/kentron)**: the kernel, [on crates.io](https://crates.io/crates/kentron). Pure
  Rust, with two optional dependencies: `rayon`, for running independent steps in parallel,
  and `ttf-parser`, for text (the `text` feature).
- **[`crates/pykentron`](crates/pykentron)**: Python bindings with a CadQuery-compatible API,
  [on PyPI](https://pypi.org/project/pykentron/). CadQuery scripts run on kentron by
  importing `pykentron` in place of `cadquery`, or unchanged with `pykentron-run`.
- **[`crates/occt-tests`](crates/occt-tests)**: OCCT's own tests, converted to run against kentron.

## Aims and scope

- **OCCT's functionality, not its API.** Names, conventions and data structures are
  kentron's own. All topology lives in one arena (`topology::Model`) addressed by typed ids,
  with no reference counting. [`docs/occt-mapping.md`](docs/occt-mapping.md) maps OCCT concepts to kentron's.
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
So far (2026-09-25 to 2026-10-29): about 252 hours of the agent's active work, 16.9 million
tokens written by it, and 9.7 billion tokens processed in all (nearly all of them its
conversation's context, read back from the cache).
<!-- /effort -->

## Tested against

- **OCCT's DRAW test scripts:** 16,311 of its 16,420 scripts in the modelling groups (the
  rest exercise viewers and data exchange). 10,245 are converted by [`tools/draw2rs.py`](tools/draw2rs.py) into
  Rust tests and YAML cases; the 6,066 that use Tcl (conditions, loops, procedures) run as
  they are, on a Tcl interpreter. Many need OCCT's external test dataset (not in this
  repo).
- **OCCT's GTests:** the computational-geometry and topology toolkits, translated by hand.
- **BRL-CAD:** NIST's MBE PMI test parts and other STEP files from its source tree, with
  OCCT's reading of them as the reference, plus its `libbg` geometry tests.
- **CadQuery's own test suite** (v2.8.0), run unchanged against pykentron: 610 of 727 tests
  pass. Most of the rest exercise OCCT's document formats (OCAF/XDE), exporters kentron
  lacks (VRML, VTK), or OCCT-specific behaviour.
- **44 CadQuery models from a private project** build on pykentron; for the 42 that return
  their parts, every part is valid and within 0.1% of CadQuery's volume.
- **748 CadQuery scripts from 113 public projects** (3 more open render windows and are not
  run): of the 739 whose results don't depend on the order a kernel lists faces in, pykentron
  builds 664 as CadQuery does, and 16 more differ only where CadQuery's result is wrong
  (2026-10-30). Timed on 2026-10-30, on the 664 it then built alike, it took 185 s in all
  against CadQuery's 1,006 s, and was faster on 633.
  The scripts, what they found and where they stand:
  [`crates/pykentron/tests/found`](crates/pykentron/tests/found/README.md).

About 12,800 tests run (with OCCT's dataset). About 1,500 of them are known gaps, marked as
expected failures ([`crates/occt-tests/expected_failures.txt`](crates/occt-tests/expected_failures.txt) for OCCT's scripts): they pass
while they fail, so `cargo test` is green, and fail once they pass, so a fix is noticed
as surely as a regression.

## Performance

A complex five-part CadQuery model, built on pykentron under each of its profiles and on
CadQuery 2.8 (OCCT), median of 5 runs:

| | wall time | CPU time | peak memory |
|---|---|---|---|
| pykentron, `fast` profile (default) | 3.3 s | 3.9 s | 406 MB |
| pykentron, `low-memory` profile | 3.9 s | 4.5 s | 305 MB |
| CadQuery 2.8 / OCCT | 37.8 s | 119.6 s | 1,089 MB |

With the default `fast` profile, that is 11× faster than CadQuery, with 31× less CPU time and 2.7×
less memory. A second complex model from the same project builds in 3.4 s (4.0 s of CPU time,
469 MB) on pykentron's `fast` profile, against 37.6 s (119.4 s, 1,103 MB) on CadQuery.

pykentron frees shapes that Python no longer refers to, with the data cached for them. Which
allocator it uses is chosen with `PYKENTRON_PROFILE`:

- `low-memory`: the system's.
- `fast` (default): mimalloc.

Rust programs using kentron directly make the same choice by declaring their allocator.

## Installing

- Rust: `cargo add kentron`.
- Python 3.10 or later: `pip install pykentron`, or `pip install "pykentron[cadquery]"` to
  run CadQuery scripts (it adds what CadQuery's own Python layers need). Wheels are built
  for Linux on x86-64-v2 (SSE4.2: Intel from 2008, AMD from 2011); elsewhere pip builds from
  the source distribution, which needs Rust.

```python
import pykentron as cq
box = cq.Workplane().box(10, 20, 30).edges("|Z").fillet(2)
```

An existing CadQuery script runs unchanged with `pykentron-run script.py`, or after
`import pykentron.shim` (before anything imports `cadquery`).

## Building and testing

- Rust 1.88 or later. The workspace builds for x86-64-v2 ([`.cargo/config.toml`](.cargo/config.toml)); on an
  older x86-64 CPU, build with an empty `RUSTFLAGS=`.
- `cargo test` runs everything. Most of OCCT's converted tests read its test dataset, an
  external download: see [`crates/occt-tests/README.md`](crates/occt-tests/README.md). Without it those tests fail.
- [`tools/run_tests.py`](tools/run_tests.py) runs the same tests with every test binary at once (about 33 s on 16
  threads, against about 49 s for `cargo test --release`, which runs them one after
  another; the slowest tests expected to fail take most of it). `--quick` leaves out the
  tests expected to fail (about 12 s); a fix that makes one pass
  then goes unnoticed.
- pykentron from source: `python3 crates/pykentron/build.py` (see
  [`crates/pykentron/README.md`](crates/pykentron/README.md)).

## More

- [`crates/pykentron/README.md`](crates/pykentron/README.md): building and using the Python bindings.
- [`crates/occt-tests/README.md`](crates/occt-tests/README.md): running OCCT's tests, and getting their data.
- [`crates/pykentron/tests/found/README.md`](crates/pykentron/tests/found/README.md): the CadQuery scripts from public projects pykentron
  is tested against, and where they stand.

## Licence

kentron and pykentron are dual-licensed under the [MIT](LICENSE-MIT) and
[Apache 2.0](LICENSE-APACHE) licences, at your option. These parts keep the licences of the
code they come from:
- [`crates/occt-tests`](crates/occt-tests), translated from OCCT's and BRL-CAD's tests: LGPL-2.1 (OCCT's with the
  Open CASCADE exception). Nothing else depends on it.
- The CadQuery code copied into pykentron
  ([`crates/pykentron/python/pykentron/cadquery_shim/cadquery`](crates/pykentron/python/pykentron/cadquery_shim/cadquery) and
  [`crates/pykentron/python/pykentron/occ_impl/exporters`](crates/pykentron/python/pykentron/occ_impl/exporters), each with its
  `LICENSE.CadQuery` and a list of the files copied in `VENDORED.txt`): Apache-2.0.
- The CadQuery scripts in [`crates/pykentron/tests/found`](crates/pykentron/tests/found), collected from public projects
  on GitHub to test pykentron against CadQuery, and the scripts derived from them in
  [`found/derived`](crates/pykentron/tests/found/derived): each keeps its original project's permissive licence (MIT, Apache-2.0,
  BSD-3-Clause, the Unlicense or CC0-1.0), copied beside the project's scripts
  (`found/<project>/LICENSE`) or, for a derived script, in its header. Scripts under other
  licences, or none, are not included;
  [`found/cases.tsv`](crates/pykentron/tests/found/cases.tsv) and [`found/passing.tsv`](crates/pykentron/tests/found/passing.tsv) give every script's source and licence.
