#!/usr/bin/env python3
"""Runs a function of a CadQuery project's module (by default `build()`) on pykentron, through
the `cadquery` shim, or on the real CadQuery (`--oracle`), and reports how far it gets: the
time taken, and for each part its volume, faces and validity. Nothing is written into the
project.

`build()` may return a shape, a dict of named shapes, or a tuple whose first item is one of
these. Run it with a Python that has the libraries the project imports (and CadQuery, for
`--oracle`):

    PYTHON crates/pykentron/scripts/cadquery_build.py --project DIR MODULE [--call build] [--oracle]

`--project` is the directory the project's modules are imported from.
"""

import argparse
import importlib
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYTHON = HERE.parent / "python"

ap = argparse.ArgumentParser()
ap.add_argument("module")
ap.add_argument("--project", required=True, help="the directory the project's modules are imported from")
ap.add_argument("--oracle", action="store_true", help="run on the real CadQuery instead, for comparison")
ap.add_argument("--call", default="build", help="the function to run (dotted: m20.m14.build)")
args = ap.parse_args()
if not args.oracle:
    sys.path[:0] = [str(PYTHON / "pykentron" / "cadquery_shim"), str(PYTHON)]

sys.path.insert(0, str(Path(args.project).expanduser().resolve()))
start = time.time()
module = importlib.import_module(args.module)
target = module
for part in args.call.split("."):
    target = getattr(target, part)
try:
    result = target()
except Exception:
    traceback.print_exc()
    print(f"FAILED after {time.time() - start:.1f} s")
    sys.exit(1)
print(f"done in {time.time() - start:.1f} s")
parts = result[0] if isinstance(result, tuple) else result
if hasattr(parts, "Volume"):
    parts = {args.module: parts}
if isinstance(parts, dict):
    for name, shape in parts.items():
        if hasattr(shape, "Volume"):
            print(f"  {name}: volume {shape.Volume():.1f}, faces {len(shape.Faces())}, valid {shape.isValid()}")
