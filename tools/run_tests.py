#!/usr/bin/env python3
"""Builds and runs the workspace's tests: every test binary of kentron (its library tests)
and of occt-tests at once, rather than one after another as `cargo test` runs them (their
slowest tests would otherwise leave most cores idle), with any other commands given
alongside.

Usage: tools/run_tests.py [--quick] [--only LIST] [--log FILE] [--also COMMAND]... [--before COMMAND]...

- `--quick` leaves out the OCCT tests expected to fail (`KENTRON_SKIP_EXPECTED_FAILURES`):
  about a third of the work, but a test starting to pass goes unnoticed.

- `--only LIST` runs only the tests a file names, one a line, each as `<binary's source>::<test>`
  (`unittests src/lib.rs::geometry::fit::tests::...`, `tests/cases.rs::occt::...`): a
  representative subset for routine runs (the fast gate's); binaries with none of them are
  not run.

- The binaries' outputs are joined into one log (default target/test-run.log) in the order
  and with the `Running <source> (<binary>)` headers that `cargo test` prints.
- Each YAML case's time is kept in target/case-times.txt for the next run, which starts
  the slowest cases first (`KENTRON_CASE_ORDER`).
- `--before` commands run after the build and before the tests (building pykentron, which
  the `--also` commands may use); `--also` commands run alongside the tests, their output
  printed after them.
- Exits with 1 if any test binary or command failed.
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TIMES = ROOT / "target" / "case-times.txt"

ap = argparse.ArgumentParser()
ap.add_argument("--quick", action="store_true", help="leave out the tests expected to fail")
ap.add_argument("--only", help="a file of the tests to run (<binary's source>::<test> a line)")
ap.add_argument("--log", default=str(ROOT / "target" / "test-run.log"))
ap.add_argument("--also", action="append", default=[], help="a shell command run alongside the tests")
ap.add_argument("--before", action="append", default=[], help="a shell command run after the build, before the tests")
args = ap.parse_args()


def timed(label, started):
    print(f"TIME {time.time() - started:.0f}s: {label}", flush=True)


# The test binaries, as cargo builds them (with their sources, for the log's headers).
started = time.time()
binaries = []
for selection in (["-p", "kentron", "--lib"], ["-p", "occt-tests"]):
    build = subprocess.run(["cargo", "test", "--release", *selection, "--no-run", "--message-format=json"], cwd=ROOT, capture_output=True, text=True)
    if build.returncode != 0:
        sys.stderr.write(build.stderr)
        sys.exit(1)
    for line in build.stdout.splitlines():
        message = json.loads(line)
        if message.get("reason") == "compiler-artifact" and message.get("executable") and message["profile"]["test"]:
            target = message["target"]
            source = str(Path(target["src_path"]).relative_to(Path(message["manifest_path"]).parent))
            # (Named as cargo names them: a library's tests as "unittests <source>".)
            label = f"unittests {source}" if "lib" in target["kind"] else source
            binaries.append((label, message["executable"], Path(message["manifest_path"]).parent))
timed("build", started)

for command in args.before:
    started = time.time()
    if subprocess.run(command, shell=True, cwd=ROOT).returncode != 0:
        sys.exit(1)
    timed(command, started)

env = dict(os.environ, KENTRON_CASE_TIMES="1", KENTRON_CASE_ORDER=str(TIMES))
if args.quick:
    env["KENTRON_SKIP_EXPECTED_FAILURES"] = "1"


# (With --only, each binary's tests from the list: the custom harnesses, whose runner takes a
# single filter, read them from a file named by KENTRON_TEST_LIST; the others take them as
# exact filters.)
only = None
if args.only:
    only = {}
    for line in Path(args.only).read_text().splitlines():
        for source, *_ in binaries:
            if line.startswith(source + "::"):
                only.setdefault(source, []).append(line[len(source) + 2:])
                break
    binaries = [b for b in binaries if b[0] in only]
CUSTOM = ("tests/cases.rs", "tests/draw_tcl.rs")


def run_binary(entry):
    source, executable, cwd = entry
    command, run_env = [executable], env
    if only is not None:
        if source in CUSTOM:
            listed = ROOT / "target" / f"only-{Path(source).stem}.txt"
            listed.write_text("".join(n + "\n" for n in only[source]))
            run_env = dict(env, KENTRON_TEST_LIST=str(listed))
        else:
            command += ["--exact", *only[source]]
    done = subprocess.run(command, cwd=cwd, env=run_env, capture_output=True, text=True)
    return entry, done.returncode, done.stdout + done.stderr


def run_command(command):
    done = subprocess.run(command, shell=True, cwd=ROOT, capture_output=True, text=True)
    return command, done.returncode, done.stdout + done.stderr


started = time.time()
with ThreadPoolExecutor(len(binaries) + len(args.also)) as pool:
    tests = [pool.submit(run_binary, b) for b in binaries]
    others = [pool.submit(run_command, c) for c in args.also]
    results = [t.result() for t in tests]
    timed("tests", started)
    commands = [o.result() for o in others]
timed("tests and commands", started)

failed = False
times = {}
with open(args.log, "w") as log:
    for (source, executable, _), code, output in results:
        log.write(f"     Running {source} ({executable})\n{output}")
        failed |= code != 0
        for m in re.finditer(r"^case time ([\d.]+)s (\S+)$", output, re.M):
            times[m.group(2)] = float(m.group(1))
summary = [(source, re.findall(r"^test result: .*$", output, re.M)) for (source, _, _), _, output in results]
passed = sum(int(m) for _, lines in summary for line in lines for m in re.findall(r"(\d+) passed", line))
failing = sum(int(m) for _, lines in summary for line in lines for m in re.findall(r"(\d+) failed", line))
print(f"tests: {passed} passed, {failing} failed (log: {args.log})")
if times:
    # (Cases not run this time, as with --quick, keep their times from before.)
    if TIMES.exists():
        for line in TIMES.read_text().splitlines():
            t, _, name = line.partition(" ")
            if name and name not in times:
                times[name] = float(t)
    TIMES.write_text("".join(f"{t:.3f} {name}\n" for name, t in sorted(times.items(), key=lambda x: -x[1])))
for command, code, output in commands:
    print(f"--- {command} (exit {code})\n{output.rstrip()}")
    failed |= code != 0
sys.exit(1 if failed else 0)
