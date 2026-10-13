#!/usr/bin/env python3
"""Copies CadQuery's pure-Python layers (no OCCT calls) into the `cadquery` shim, unchanged.

CadQuery's `Workplane`, `Sketch`, selectors, `cqgi`, free functions and assemblies are written
on its `occ_impl` layer, which the shim provides from pykentron; this copies those modules from
a CadQuery release (Apache-2.0: its licence is copied with them). Run it again to update them
after changing the release.

Usage:
    python3 crates/pykentron/scripts/vendor_cadquery.py [--source DIR]
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHIM = HERE.parent / "python" / "cadquery_shim" / "cadquery"
EXPORTERS = HERE.parent / "python" / "pykentron" / "occ_impl" / "exporters"
REPOSITORY = "https://github.com/CadQuery/cadquery.git"
RELEASE = "v2.8.0"

# The modules copied (relative to CadQuery's `cadquery` package): those with no direct OCCT
# calls; `__init__.py` and the rest of `occ_impl` are the shim's own.
VENDORED = [
    "assembly.py",
    "cq.py",
    "cqgi.py",
    "func.py",
    "hull.py",
    "selectors.py",
    "sketch.py",
    "types.py",
    "units.py",
    "utils.py",
    "plugins/__init__.py",
    "occ_impl/shape_protocols.py",
    "occ_impl/sketch_solver.py",
    "occ_impl/solver.py",
]

# CadQuery's mesh writers (pure Python, on shapes' tessellations), into pykentron's
# `occ_impl.exporters`, where its other exporters are.
VENDORED_EXPORTERS = ["amf.py", "json.py", "threemf.py"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", help="a CadQuery checkout at the release (fetched if not given)")
    args = ap.parse_args()
    if args.source:
        src = Path(args.source)
    else:
        src = HERE.parent.parent.parent / "target" / "cadquery-tests" / f"cadquery-{RELEASE}"
        if not (src / "cadquery").is_dir():
            subprocess.run(["git", "clone", "--quiet", "--depth", "1", "--branch", RELEASE, REPOSITORY, str(src)], check=True)
    for rel in VENDORED:
        target = SHIM / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src / "cadquery" / rel, target)
    for rel in VENDORED_EXPORTERS:
        shutil.copyfile(src / "cadquery" / "occ_impl" / "exporters" / rel, EXPORTERS / rel)
    for target, files, prefix in ((SHIM, VENDORED, ""), (EXPORTERS, VENDORED_EXPORTERS, "occ_impl/exporters/")):
        shutil.copyfile(src / "LICENSE", target / "LICENSE.CadQuery")
        (target / "VENDORED.txt").write_text(
            f"Copied unchanged from CadQuery {RELEASE} ({REPOSITORY}), Apache License 2.0 "
            "(LICENSE.CadQuery), by crates/pykentron/scripts/vendor_cadquery.py:\n" + "".join(f"  {prefix}{r}\n" for r in files)
        )
    print(f"copied {len(VENDORED) + len(VENDORED_EXPORTERS)} modules from CadQuery {RELEASE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
