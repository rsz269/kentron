#!/usr/bin/env python3
"""Builds every module of a CadQuery project on pykentron and compares the parts with the real
CadQuery's (the oracle): each part valid, its volume within 0.1%. Builds run in parallel,
each in a process of its own, through scripts/cadquery_build.py.

    crates/pykentron/scripts/cadquery_sweep.py --project DIR --python PYTHON [--cwd DIR]
        [--exclude MODULE,...] [--oracle-cache FILE] [MODULE ...]

- `--project`: the directory the project's modules are imported from; the modules swept are
  those given, or else every one there defining `build()` at its top level.
- `--python`: a Python with the project's libraries and CadQuery (default: this one).
- `--cwd`: the directory builds run in (some projects read files relative to it).
- `--oracle-cache`: a JSON file keeping the oracle's volumes (they do not change): only
  modules missing from it are built on CadQuery.

Prints each mismatch, then a summary; exits with 1 if any module fails to build or differs.
"""

import argparse
import json
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BUILD = Path(__file__).resolve().parent / "cadquery_build.py"

ap = argparse.ArgumentParser()
ap.add_argument("modules", nargs="*")
ap.add_argument("--project", required=True)
ap.add_argument("--python", default=sys.executable)
ap.add_argument("--cwd", default=".")
ap.add_argument("--exclude", default="", help="modules to leave out, comma-separated")
ap.add_argument("--oracle-cache", type=Path)
ap.add_argument("--jobs", type=int, default=8)
args = ap.parse_args()

project = Path(args.project).expanduser().resolve()
excluded = set(filter(None, args.exclude.split(",")))
modules = args.modules or sorted(
    p.stem for p in project.glob("*.py") if re.search(r"^def build\(", p.read_text(), re.M)
)
modules = [m for m in modules if m not in excluded]


def run(module: str, oracle: bool) -> tuple[str, bool, dict]:
    command = [args.python, str(BUILD), "--project", str(project), module] + (["--oracle"] if oracle else [])
    out = subprocess.run(command, cwd=args.cwd, capture_output=True, text=True, timeout=900).stdout
    parts = {m[0]: (float(m[1]), m[2] == "True") for m in re.findall(r"^\s+(\S+): volume (\S+), faces \d+, valid (\w+)", out, re.M)}
    return module, "done in" in out, parts


cache = args.oracle_cache
oracle = json.loads(cache.read_text()) if cache and cache.exists() else {}
with ThreadPoolExecutor(args.jobs) as ex:
    for module, ok, parts in ex.map(lambda m: run(m, True), [m for m in modules if m not in oracle]):
        if ok:
            oracle[module] = {name: volume for name, (volume, _) in parts.items()}
if cache:
    cache.write_text(json.dumps(oracle, indent=1, sort_keys=True) + "\n")

built = matched = 0
with ThreadPoolExecutor(args.jobs) as ex:
    for module, ok, parts in ex.map(lambda m: run(m, False), modules):
        built += ok
        want = oracle.get(module)
        good = ok and want is not None and set(parts) == set(want) and all(
            parts[k][1] and abs(parts[k][0] - want[k]) <= 1e-3 * max(abs(want[k]), 1.0) for k in want
        )
        matched += good
        if not good:
            print(f"MISMATCH {module}: built={ok} parts={parts} oracle={want}")
print(f"{built}/{len(modules)} build, {matched} match the oracle (volumes within 0.1%, all valid)")
sys.exit(0 if matched == len(modules) else 1)
