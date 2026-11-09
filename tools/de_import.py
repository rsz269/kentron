#!/usr/bin/env python3
"""Copies OCCT's data-exchange tests (tests/de: STEP and IGES files read, written and read
again, their statistics compared with reference data) into crates/occt-tests/de/, for the
draw_de harness (crates/occt-tests/tests/draw_de.rs).

Only the cases whose file is in OCCT's public dataset (or OCCT's source tree) are copied:
the others read models OCCT keeps private and can never run. They are excluded for good,
listed in crates/occt-tests/de/excluded.txt with the file each needs.

    tools/de_import.py [--occt ~/OCCT] [--data DIR ...]

The data is looked for as the harness looks for it (find_data_file): by exact name under
$OCCT_TEST_DATA, ~/occt-data, $OCCT_ROOT/data and ~/OCCT/data.
"""

import argparse
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "crates", "occt-tests", "de")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--occt", default=os.path.expanduser("~/OCCT"))
    ap.add_argument("--data", action="append", default=[])
    args = ap.parse_args()
    dirs = list(args.data) + [d for d in os.environ.get("OCCT_TEST_DATA", "").split(os.pathsep) if d]
    dirs += [os.path.expanduser("~/occt-data")]
    if os.environ.get("OCCT_ROOT"):
        dirs.append(os.path.join(os.environ["OCCT_ROOT"], "data"))
    dirs.append(os.path.join(args.occt, "data"))
    have = set()
    for d in dirs:
        for _, _, files in os.walk(d):
            have.update(files)
    source = os.path.join(args.occt, "tests", "de")
    try:
        commit = subprocess.run(["git", "-C", args.occt, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = ""
    if os.path.isdir(OUT):
        for entry in os.listdir(OUT):
            path = os.path.join(OUT, entry)
            if os.path.isdir(path):
                shutil.rmtree(path)
    os.makedirs(OUT, exist_ok=True)
    copied, excluded = 0, []
    for grid in sorted(os.listdir(source)):
        gdir = os.path.join(source, grid)
        if not os.path.isdir(gdir):
            continue
        for case in sorted(os.listdir(gdir)):
            if case in ("begin", "end"):
                continue
            text = open(os.path.join(gdir, case), errors="replace").read()
            m = re.search(r"set filename\s+(\S+)", text)
            name = m.group(1) if m else ""
            if name not in have:
                excluded.append((grid, case, name))
                continue
            os.makedirs(os.path.join(OUT, grid), exist_ok=True)
            shutil.copy(os.path.join(gdir, case), os.path.join(OUT, grid, case))
            copied += 1
    with open(os.path.join(OUT, "excluded.txt"), "w") as f:
        f.write("# OCCT's tests/de cases left out for good: each reads a file that is not in OCCT's public\n")
        f.write("# dataset (OCCT keeps most of its test models private), so it can never run.\n")
        f.write(f"# Written by tools/de_import.py from OCCT {commit}.\n")
        f.write("# grid/case\tfile\n")
        for grid, case, name in excluded:
            f.write(f"{grid}/{case}\t{name}\n")
    print(f"copied {copied} cases, excluded {len(excluded)} (no data)")


if __name__ == "__main__":
    main()
