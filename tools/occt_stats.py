#!/usr/bin/env python3
"""How OCCT's DRAW scripts stand in kentron's suite, for the README's figures.

Each script of the modelling groups is counted once, wherever it runs (a Rust test or YAML
case converted by draw2rs, or the script itself on the Tcl interpreter), as one of:

- no data: it reads a file that is not in OCCT's public dataset (nor OCCT's source tree):
  OCCT keeps most of its test models private, so these never run, and count neither as
  passing nor as failing;
- known OCCT failure: OCCT's own script declares it fails on this platform (a TODO), and it
  is not run;
- expected failure: it runs and fails (listed in expected_failures.txt, or a YAML case
  `expect: not_done`): a known gap;
- passing: it runs and passes.

The data files are looked for as the harness does (crates/occt-tests/src/draw/mod.rs,
find_data_file): by exact name under $OCCT_TEST_DATA, ~/occt-data, $OCCT_ROOT/data,
~/OCCT/data and BRL-CAD's source tree. Run with the public dataset unpacked in ~/occt-data
(as CI has it).
"""

import glob
import os
import re
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "crates", "occt-tests")


def data_index():
    dirs = [d for d in os.environ.get("OCCT_TEST_DATA", "").split(os.pathsep) if d]
    home = os.path.expanduser("~")
    dirs.append(os.path.join(home, "occt-data"))
    if os.environ.get("OCCT_ROOT"):
        dirs.append(os.path.join(os.environ["OCCT_ROOT"], "data"))
    dirs.append(os.path.join(home, "OCCT", "data"))
    brlcad = os.environ.get("BRLCAD_SRC", os.path.join(home, "brlcad"))
    dirs += [os.path.join(brlcad, "db"), os.path.join(brlcad, "src", "libbrep", "tests")]
    names = set()
    for d in dirs:
        for _, _, files in os.walk(d):
            names.update(files)
    return names


TODO = re.compile(r'^\s*puts\s+"(TODO\s+\S+\s+([^:"]+):[^"]*)"')


def known_failure(text):
    """As occt_tests::draw::tcl::known_failure: a TODO for this platform outside braces."""
    depth = 0
    for line in text.splitlines():
        if depth == 0:
            m = TODO.match(line)
            if m and any(p in ("ALL", "All", "Linux", "Debian", "Ubuntu") for p in re.split(r"[ ,]", m.group(2))):
                return True
        depth += line.count("{") - line.count("}")
    return False


def main():
    have = data_index()
    expected = {l.strip() for l in open(os.path.join(ROOT, "expected_failures.txt")) if l.startswith("tests/")}
    status = {}  # script ("tests/group/grid/case") -> status

    # The scripts run on the Tcl interpreter.
    data_file = re.compile(r"locate_data_file\s+([^\]\s}]+)")
    for path in glob.glob(os.path.join(ROOT, "scripts", "*", "*", "*")):
        group, grid, case = path.split(os.sep)[-3:]
        if case in ("begin", "end"):
            continue
        name = f"tests/{group}/{grid}/{case}"
        text = open(path, errors="replace").read()
        missing = [f for f in data_file.findall(text) if "$" not in f and f not in have]
        if missing:
            status[name] = "no data"
        elif known_failure(text):
            status[name] = "known OCCT failure"
        elif name in expected:
            status[name] = "expected failure"
        else:
            status[name] = "passing"

    # The converted Rust tests.
    test = re.compile(r"/// `(tests/[^`]+)`\n(?:///[^\n]*\n)*((?:#\[[^\n]*\]\n)*)fn ")
    for path in glob.glob(os.path.join(ROOT, "tests", "draw_*", "*.rs")):
        text = open(path).read()
        for m in test.finditer(text):
            name, attributes = m.group(1), m.group(2)
            if "missing OCCT test data" in attributes:
                status[name] = "no data"
            elif "known failure in OCCT" in attributes:
                status[name] = "known OCCT failure"
            elif name in expected or "should_panic" in attributes:
                status[name] = "expected failure"
            else:
                status[name] = "passing"

    # The converted YAML cases (cases/occt/draw), as occt_tests::cases::load reads them.
    try:
        import yaml
    except ImportError:
        sys.exit("needs PyYAML")
    for path in glob.glob(os.path.join(ROOT, "cases", "occt", "draw", "**", "*.yaml"), recursive=True):
        for doc in yaml.safe_load_all(open(path)):
            name = (doc or {}).get("ported_from", {}).get("file")
            if not name:
                continue
            loads = []
            for step in doc.get("build") or []:
                if isinstance(step, dict):
                    for value in step.values():
                        load = value.get("load") if isinstance(value, dict) else None
                        if isinstance(load, dict):
                            loads += [load[k] for k in ("brep", "step") if isinstance(load.get(k), str)]
            if doc.get("skip"):
                status[name] = "known OCCT failure"
            elif any(f not in have for f in loads):
                status[name] = "no data"
            elif doc.get("expect") == "not_done":
                status[name] = "expected failure"
            else:
                status[name] = "passing"

    counts = {}
    for s in status.values():
        counts[s] = counts.get(s, 0) + 1
    run = counts.get("passing", 0) + counts.get("expected failure", 0)
    print(f"scripts: {len(status)}")
    for s in ("no data", "known OCCT failure", "expected failure", "passing"):
        print(f"  {s}: {counts.get(s, 0)}")
    print(f"run: {run}, passing {counts.get('passing', 0)} ({100 * counts.get('passing', 0) / max(run, 1):.1f}%)")
    if "--by-group" in sys.argv:
        groups = {}
        for name, s in status.items():
            g = groups.setdefault(name.split("/")[1], {})
            g[s] = g.get(s, 0) + 1
        for g, c in sorted(groups.items(), key=lambda x: -sum(x[1].values())):
            r = c.get("passing", 0) + c.get("expected failure", 0)
            print(f"  {g:12} run {r:5}  passing {c.get('passing', 0):5}  expected failures {c.get('expected failure', 0):4}  no data {c.get('no data', 0):5}  OCCT TODO {c.get('known OCCT failure', 0):4}")


if __name__ == "__main__":
    main()
