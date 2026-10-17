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
