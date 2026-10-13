#!/usr/bin/env python3
"""Builds the native module and places it in the `pykentron` package (no maturin needed).

Usage: crates/pykentron/build.py [--python PATH] [--debug | --portable | --quick]
Optimized builds use the `dist` profile (link-time optimization across crates) and AVX2/FMA
(x86-64-v3, CPUs from 2013 on); `--portable` builds for any x86-64; `--quick` builds with the
`release` profile the tests use (reusing their build of kentron: for checking changes, not
for shipping). The package is then
importable with crates/pykentron/python on PYTHONPATH (and the
`cadquery` shim with crates/pykentron/python/cadquery_shim too).
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent

ap = argparse.ArgumentParser()
ap.add_argument("--python", default=sys.executable, help="the Python the module is built for")
ap.add_argument("--debug", action="store_true")
ap.add_argument("--portable", action="store_true", help="no instructions beyond baseline x86-64")
ap.add_argument("--quick", action="store_true", help="the tests' `release` profile, reusing their build of kentron")
args = ap.parse_args()

env = dict(os.environ, PYO3_PYTHON=args.python)
profile = "debug" if args.debug else "release" if args.quick else "dist"
# (The workspace's .cargo/config.toml builds for x86-64-v3 already; a portable build overrides it.)
if args.portable and "RUSTFLAGS" not in os.environ:
    env["RUSTFLAGS"] = "-C target-cpu=x86-64"
subprocess.run(["cargo", "build", "-p", "pykentron", "--profile", "dev" if args.debug else profile], cwd=ROOT, env=env, check=True)
# (Where cargo put it: CARGO_TARGET_DIR, relative to the workspace, else its target/.)
built = ROOT / os.environ.get("CARGO_TARGET_DIR", "target") / profile / "lib_kentron.so"
target = HERE / "python" / "pykentron" / "_kentron.abi3.so"
# Replaced by a rename, so that running interpreters keep the library they loaded.
staging = target.with_name(target.name + ".new")
shutil.copyfile(built, staging)
os.replace(staging, target)
print(f"built {target}")
