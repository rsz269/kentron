#!/usr/bin/env python3
"""Runs CadQuery's own test suite against pykentron, unchanged.

The tests are fetched from CadQuery's repository at a pinned release, into a cache, and run
with pytest in a virtual environment of their own (made with `uv`), with the `cadquery` shim
first on the path, so that `import cadquery` loads pykentron. Each test's outcome is compared
with the baseline (`crates/pykentron/tests/cadquery_baseline.txt`, the tests that fail): the run
fails on a test failing that is not in it, and reports those of it that pass.

Usage:
    python3 crates/pykentron/scripts/cadquery_tests.py [--update-baseline] [pytest arguments...]

Build the extension first (`python3 crates/pykentron/build.py`).
"""

import argparse
import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYKENTRON = HERE.parent
ROOT = PYKENTRON.parent.parent
BASELINE = PYKENTRON / "tests" / "cadquery_baseline.txt"

# The CadQuery release whose tests run, and the Python packages they import besides
# cadquery and OCP (which the shim provides).
REPOSITORY = "https://github.com/CadQuery/cadquery.git"
RELEASE = "v2.8.0"
PYTHON = "3.12"
PACKAGES = [
    "pytest", "pytest-timeout", "pytest-xdist", "numpy", "scipy", "docutils", "ezdxf", "path",
    # (Those of the CadQuery layers the shim runs.)
    "pyparsing", "multimethod", "runtype", "typing_extensions", "nlopt", "casadi",
]


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, **kw)


def sources(cache: Path) -> Path:
    """CadQuery's sources at the pinned release, fetched once."""
    src = cache / f"cadquery-{RELEASE}"
    if not (src / "tests").is_dir():
        run(["git", "clone", "--quiet", "--depth", "1", "--branch", RELEASE, REPOSITORY, str(src)])
    return src


def environment(cache: Path) -> Path:
    """A virtual environment with the tests' packages, made once; its Python."""
    venv = cache / "venv"
    python = venv / "bin" / "python"
    stamp = venv / "packages.txt"
    wanted = "\n".join([PYTHON] + PACKAGES)
    if not python.exists() or not stamp.exists() or stamp.read_text() != wanted:
        shutil.rmtree(venv, ignore_errors=True)
        run(["uv", "venv", "--quiet", "--python", PYTHON, str(venv)])
        run(["uv", "pip", "install", "--quiet", "--python", str(python), *PACKAGES])
        stamp.write_text(wanted)
    return python


def outcomes(junit: Path) -> dict:
    """Test id -> passed / failed / skipped, from pytest's JUnit report."""
    out = {}
    for case in ET.parse(junit).getroot().iter("testcase"):
        name = f"{case.get('classname')}::{case.get('name')}"
        if case.find("failure") is not None or case.find("error") is not None:
            out[name] = "failed"
        elif case.find("skipped") is not None:
            out[name] = "skipped"
        else:
            out[name] = "passed"
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", default=str(ROOT / "target" / "cadquery-tests"), help="where sources and the environment are kept")
    ap.add_argument("--update-baseline", action="store_true", help="record this run's failures as the baseline")
    ap.add_argument("--timeout", type=int, default=120, help="seconds per test")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2), help="tests run at once (pytest-xdist workers; 1 runs them in turn)")
    args, pytest_args = ap.parse_known_args()

    cache = Path(args.cache)
    cache.mkdir(parents=True, exist_ok=True)
    src = sources(cache)
    python = environment(cache)

    # The tests alone, apart from CadQuery's own package (which would be imported instead of
    # the shim from beside them).
    suite = cache / "suite"
    shutil.rmtree(suite, ignore_errors=True)
    shutil.copytree(src / "tests", suite / "tests")
    for f in ("pytest.ini", "conftest.py"):
        if (src / f).exists():
            shutil.copy(src / f, suite / f)

    junit = cache / "junit.xml"
    junit.unlink(missing_ok=True)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(PYKENTRON / "python" / "pykentron" / "cadquery_shim"), str(PYKENTRON / "python")])
    cmd = [str(python), "-m", "pytest", "tests", "-q", "-p", "no:cacheprovider", "--continue-on-collection-errors", f"--timeout={args.timeout}", f"--junitxml={junit}", "-o", "junit_family=xunit1", *(["-n", str(args.jobs)] if args.jobs > 1 else []), *pytest_args]
    print("+", " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, cwd=suite, env=env)
    if not junit.exists():
        print("pytest wrote no report (the suite did not load)", file=sys.stderr)
        return 2

    results = outcomes(junit)
    failed = sorted(t for t, o in results.items() if o == "failed")
    counts = {o: sum(1 for v in results.values() if v == o) for o in ("passed", "failed", "skipped")}
    print(f"\n{len(results)} tests: {counts['passed']} passed, {counts['failed']} failed, {counts['skipped']} skipped")

    if args.update_baseline:
        BASELINE.write_text(f"# CadQuery {RELEASE} tests that fail on pykentron (cadquery_tests.py --update-baseline)\n" + "".join(t + "\n" for t in failed))
        print(f"baseline updated: {len(failed)} failing tests")
        return 0
    known = set()
    if BASELINE.exists():
        known = {l.strip() for l in BASELINE.read_text().splitlines() if l.strip() and not l.startswith("#")}
    new = [t for t in failed if t not in known]
    fixed = sorted(t for t in known if results.get(t) == "passed")
    for t in fixed:
        print(f"now passes: {t}")
    for t in new:
        print(f"NEW FAILURE: {t}")
    print(f"{len(new)} new failures, {len(fixed)} now passing (of {len(known)} in the baseline)")
    return 1 if new else 0


if __name__ == "__main__":
    sys.exit(main())
