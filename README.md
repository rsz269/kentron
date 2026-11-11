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

kentron models solids as faces bounded by edges and vertices, over exact curves and surfaces
(analytic and NURBS). It is inspired by Open CASCADE Technology (OCCT) and written
independently: OCCT is not a dependency.

- Primitives, extrusions, revolutions, lofts and pipes.
- Booleans, fillets, chamfers, drafts, offsets and thick solids.
- Filling surfaces, removing features, sewing and healing.
- Meshing, and text from TrueType and OpenType fonts.
- BRep, STEP and STL files.

| Crate | |
|---|---|
| [`kentron`](crates/kentron) | The kernel, [on crates.io](https://crates.io/crates/kentron). Pure Rust. Optional dependencies: `rayon` (parallel steps) and `ttf-parser` (text). |
| [`pykentron`](crates/pykentron) | Python bindings with a CadQuery-compatible API, [on PyPI](https://pypi.org/project/pykentron/). Its own assembly and sketch constraint solvers replace casadi and nlopt. |
| [`occt-tests`](crates/occt-tests) | OCCT's tests, converted to run against kentron. |

## Status

- Pre-1.0.
- Names, conventions and data structures are kentron's own, not OCCT's.
  [`docs/occt-mapping.md`](docs/occt-mapping.md) maps OCCT's concepts to them.
- All topology lives in one arena (`topology::Model`) addressed by typed ids, without
  reference counting.
- Out of scope for now: visualization, hidden-line removal, OCAF/XDE documents.
- Not done yet: evolved solids, fillets at saddle corners. Rare pathological inputs can run
  very slowly or fail to converge. Known gaps: [`docs/known-gaps.md`](docs/known-gaps.md).

## How it was made

kentron is 99% vibe coded: written almost entirely by an AI coding agent (Claude, in Claude
Code), with a human setting the direction and judging the results. Every change is checked
against the tests below before it is committed.

<!-- effort -->
So far (2026-09-25 to 2026-11-11): about 364 hours of the agent's active work, 25.6 million
tokens written by it, and 13.6 billion tokens processed in all (nearly all of them its
conversation's context, read back from the cache).
<!-- /effort -->

## Tests

About 13,430 tests run. About 1,150 of them are known gaps, marked as expected failures, so
`cargo test` is green.

- **OCCT's DRAW scripts:** 9,649 of 10,796 pass (89%). Scripts that need OCCT's unpublished
  test models are not run.
- **OCCT's GTests** for geometry and topology, translated by hand.
- **USD BRep's tests:** 212 cases, translated.
- **BRL-CAD:** STEP files from its source tree and its `libbg` geometry tests.
- **CadQuery 2.8.0's test suite,** unchanged: 698 of 699 pass. 68 tests of OCCT internals
  are left out ([list](crates/pykentron/tests/cadquery_excluded.txt)).
- **809 CadQuery scripts from 154 public projects:** of the 795 whose results don't depend
  on a kernel's face order, 719 build as on CadQuery, 27 differ only where CadQuery's result
  is wrong, and 49 fail ([details](crates/pykentron/tests/found/README.md)).
- **44 CadQuery models from a private project** build. The 42 that return their parts match
  CadQuery's volumes within 0.1%, every part valid.
- **Own tests:** 433 unit tests, and 58 cases compared with CadQuery.

## Performance

A complex five-part CadQuery model, median of 5 runs:

| | wall time | CPU time | peak memory |
|---|---|---|---|
| pykentron, `fast` profile (default) | 3.5 s | 3.9 s | 381 MB |
| pykentron, `low-memory` profile | 4.2 s | 4.6 s | 281 MB |
| CadQuery 2.8 / OCCT | 41.2 s | 125.4 s | 1,075 MB |

- The 719 public CadQuery scripts that build as on CadQuery (2026-11-10): 268 s in all
  against CadQuery's 1,218 s, faster on 689.
- Assembly and sketch constraints: 11–34× faster than CadQuery's solvers, with the
  same solutions.
- `PYKENTRON_PROFILE` picks the allocator: `fast` (mimalloc, the default) or `low-memory`
  (the system's).

## Installing

- Rust: `cargo add kentron`.
- Python 3.10 or later: `pip install pykentron`, or `pip install "pykentron[cadquery]"` to
  run CadQuery scripts.
- Wheels are built for Linux (x86-64-v2, ARM64), macOS (ARM64, x86-64) and Windows (x86-64).
  Elsewhere pip builds from source, which needs Rust.

```python
import pykentron as cq
box = cq.Workplane().box(10, 20, 30).edges("|Z").fillet(2)
```

An existing CadQuery script runs unchanged with `pykentron-run script.py`, or after
`import pykentron.shim` (before anything imports `cadquery`).

## Building and testing

- Rust 1.88 or later. The workspace builds for x86-64-v2
  ([`.cargo/config.toml`](.cargo/config.toml)); on an older x86-64 CPU, build with an empty
  `RUSTFLAGS=`.
- `cargo test` runs everything. Most of OCCT's converted tests need its test dataset, a
  separate download: [`crates/occt-tests/README.md`](crates/occt-tests/README.md).
- [`tools/run_tests.py`](tools/run_tests.py) runs the test binaries in parallel (about 35 s
  on 16 threads). `--quick` skips the expected failures (about 17 s).
- pykentron from source: `python3 crates/pykentron/build.py`
  ([`crates/pykentron/README.md`](crates/pykentron/README.md)).

## Licence

[AGPL-3.0-or-later](LICENSE). Third-party parts keep their own licences:

| Part | Origin | Licence |
|---|---|---|
| [`constraints/slsqp.rs`](crates/kentron/src/constraints/slsqp.rs) | translated from NLopt's SLSQP | BSD-3-Clause and MIT ([`LICENSE-SLSQP`](LICENSE-SLSQP)) |
| [`cadquery_shim/cadquery`](crates/pykentron/python/pykentron/cadquery_shim/cadquery), [`occ_impl/exporters`](crates/pykentron/python/pykentron/occ_impl/exporters) | copied from CadQuery 2.8.0 | Apache-2.0 |
| [`crates/occt-tests`](crates/occt-tests) | translated from OCCT's and BRL-CAD's tests | LGPL-2.1 (OCCT's with the Open CASCADE exception) |
| [`occt-tests/cases/usd_brep`](crates/occt-tests/cases/usd_brep) | translated from USD BRep's tests | Apache-2.0 |
| [`pykentron/tests/found`](crates/pykentron/tests/found) | CadQuery scripts from public projects | each project's own: MIT, Apache-2.0, BSD-3-Clause, Unlicense or CC0-1.0 ([`cases.tsv`](crates/pykentron/tests/found/cases.tsv)) |
