# CadQuery scripts from public projects

pykentron is tested against CadQuery 2.8 on a corpus of 501 CadQuery scripts from 113
public projects on GitHub, collected from 2026-10-18 to 2026-10-23. Each script runs on both,
and their results are compared. The scripts are a mix: when they were found, some failed on
pykentron, built differently or ran much slower, and the rest built as on CadQuery. Many of
the first kind have been fixed since, and some differences turned out not to be pykentron's
(below).

## What is here

- `<project>/`: the scripts of projects under permissive licences, copied with the project's
  licence beside them and the project's own modules they import (344 scripts from 69
  projects). Larger libraries some scripts import are not copied: cq_warehouse, cqterrain,
  cadqueryhelper, cqindustry, cqfantasy, cq_gears, cqMore, cq-electronics and LatticeQuery.
  Install them without their dependencies (`pip install --no-deps`), so that they cannot
  replace CadQuery.
- [`cases.tsv`](cases.tsv): the 211 scripts that pykentron failed on (81), built differently
  (110) or ran substantially slower (20) when they were found: project, script, source repository and
  commit, licence, kind, suspected area, both kernels' wall times, and the first differences
  or the error. It records them as found and is not updated as they are fixed.
- [`passing.tsv`](passing.tsv): the 359 scripts that built as on CadQuery when last checked
  for it (2026-10-20, and 2026-10-23 for the last round), with both kernels' wall and CPU times
  (CadQuery's from the hunt, in a fresh process each; it runs several threads, so its CPU
  time is often over its wall time). 69 scripts are in both lists: found failing or
  differing, then fixed by the time this list was made.
- [`derived/`](derived): scripts derived from corpus scripts whose results depend on a
  kernel's face order (below), with that step made explicit. They are not among the 501.
- [`harness/`](harness): the tools that run and compare the scripts.

Scripts under other licences, or none, are not copied: [`cases.tsv`](cases.tsv) and
[`passing.tsv`](passing.tsv) still list them, with the licence as `not copied (...)`.

## Where they stand

On 2026-10-27:

| | scripts |
|---|---|
| in the corpus | 501 |
| left out: the result depends on face order (below) | 6 |
| left out: CadQuery itself fails | 0 |
| compared | 495 |
| build as on CadQuery | 433 |
| differ where CadQuery's result is wrong (below) | 14 |
| differ | 39 |
| fail on pykentron | 9 |
| **pykentron's failures** | **48** |

Of the 191 scripts that failed or built differently when found, 123 now build as on
CadQuery, and 14 more differ only where CadQuery's result is wrong.

Some differences are not pykentron's:

- **CadQuery's result is wrong**, by a measure independent of both kernels or by OCCT's own
  answer to the same question asked another way: OCCT's volume of its own thread solids (the
  contrib thread, a threaded adapter: meshed, they enclose what pykentron's do); a fillet OCCT
  makes off the rolling ball (cadquery-models' lid fixer: its face is 24.05 where the exact
  surface's is 22.96); unions OCCT makes a fifth too small (the Sabrent hub stands: Monte
  Carlo sampling of the operands agrees with pykentron); a union of a compound whose own
  solids overlap or touch, which OCCT leaves as separate solids, counting the overlaps twice
  (cqterrain's buildings, helmecke's exercise 005: given the solids as arguments of their
  own, OCCT makes pykentron's result); a loft inside out (negative volume); dovetails dropped
  by `union(glue=True)`, which assumes shapes only touch; a fuzzy union coming out empty (the
  spindle vacuum shroud); a fillet that fails on OCCT inside a script's `try`, so that
  CadQuery's result lacks it.
- **The script depends on face or edge order**, which OCCT does not define and kentron does
  not copy: an asymmetric chamfer measures its first length on the edge's "first" face, a
  workplane on several faces takes the first one's normal, an edge picked by index, a loop
  over a set of shapes. Neither kernel is wrong; these are left out of the comparison, and
  where the rest of the model is worth comparing, [`derived/`](derived) has it with the
  order-dependent step made explicit (as OCCT takes it).

## How they were collected

- **2026-10-18:** the first 20 CadQuery projects in which CadQuery itself built anything.
- **2026-10-19:** 57 more from GitHub's repository search (`cadquery`, `topic:cadquery`, and
  `cadquery` with gridfinity, keyboard, enclosure, printable, parametric and the like),
  leaving out lists, AI pipelines and datasets: 382 scripts, of which CadQuery built 187.
- **2026-10-20:** 10 more (`cadquery` with 3d print, parametric, stl, design and the like): 61
  scripts, of which CadQuery built 44.
- **2026-10-23:** the 400 most-starred Python projects of a wider search (1,959 repositories
  mentioning CadQuery; lists, CadQuery itself, its editors and build123d left out), up to 30
  scripts each; CadQuery built about 530. 149 were taken, at most 15 a project, to bring the
  corpus to its size: 49 that pykentron failed on or built differently, 15 where it was
  slower or less than 3× faster, and 85 at random among those it built alike.

Left out throughout: scripts using unseeded randomness (no two runs agree), scripts calling
OCCT's Python bindings (OCP) beyond the shapes CadQuery hands them (they test OCCT's API, not
CadQuery's), scripts writing to a directory that does not exist or reading a file on their
author's disk, and scripts showing nothing to compare.

## How they are compared

[`harness/runner.py`](harness/runner.py) runs a script as `__main__` in a fresh, isolated
process, on pykentron's shim or on CadQuery, standing in for CQ-editor's `show_object` and the VS Code viewer's
module (`ocp_vscode`). It records what the script shows, exports, or leaves at module level:
volume, area, bounding box, validity, and solid and face counts. The volume is its solids'
(OCCT gives wires and faces in a compound volumes of their own); CadQuery's bounding box is
OCCT's tight one (its default is grown by tolerances and control polygons).

A script differs if pykentron raises where CadQuery does not; a volume or area differs by
more than 0.1%, or a bounding box by more than 1% of its size; the solid count differs; or
pykentron's result is invalid where CadQuery's is valid. (Exports to temporary files, named
at random, are matched by their order.) It is slower if pykentron takes
more than twice as long and 0.5 s more.

## Running them

- [`harness/hunt.py`](harness/hunt.py) `PROJECT_DIR [GLOB...]` runs up to 30 of a
  project's scripts on both kernels (120 s each) and compares them;
  [`harness/analyze.py`](harness/analyze.py) gives each flagged script a suspected area.
- [`harness/timing.py`](harness/timing.py) `--python PYTHON [FILTER...] [--repeat N]` runs
  the copied scripts of `passing.tsv` on pykentron again, each in a fresh process, and prints their times against
  CadQuery's. `PYTHON` needs the libraries above and numpy; CadQuery itself is not used.
