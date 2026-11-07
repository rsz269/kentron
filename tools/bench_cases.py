#!/usr/bin/env python3
"""Times YAML cases on kentron, OCCT and SMLib alike, inside the kernels' calls only.

Rounds alternate between the kernels (other work on the machine then falls on all); each case
runs REPEATS times a round, and each step's least time over all of them counts. Left out of
the times: parsing cases, picking sub-shapes, and the copies of operands SMLib's operations
need, consuming them as they do. Only cases all three finish count in the comparisons.

  cargo build --profile dist -p occt-tests --bin bench_cases
  tools/bench_cases.py run --kentron target/dist/bench_cases --smlib-python P --occt-python P --out OUT.json [--pin CPU]
                           [--parallel] [--rounds N] FILE...
  tools/bench_cases.py report OUT.json [--clock cpu|wall]

--pin runs each on that processor alone, one thread; --parallel lets each use its threads
(kentron's pool, OCCT's parallel booleans; SMLib has no such switch).
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import crosscheck_cases as cc

TIMED = """box cylinder cone sphere torus wedge translate rotate scale mirror union difference intersection section compound
fillet fillet_law chamfer remove_faces thick_solid offset_shape offset_wire face extrude extrude_drafted revolve pipe loft
to_nurbs helix wire vertex segment curve_edge volume area length centroid inertia bounds valid counts distance contains""".split()
CHECKS = {"volume", "area", "length", "centroid", "inertia", "bounds", "valid", "counts", "distance", "contains"}


class Timed:
    """A kernel whose operations are clocked, less what it spends on copies of operands."""

    def __init__(self, kernel):
        self.kernel, self.steps = kernel, []

    def __getattr__(self, name):
        attribute = getattr(self.kernel, name)
        if name not in TIMED:
            return attribute

        def call(*args):
            spare = self.kernel.spare = [0, 0]
            c0, w0 = time.process_time_ns(), time.perf_counter_ns()
            try:
                return attribute(*args)
            finally:
                c1, w1 = time.process_time_ns(), time.perf_counter_ns()
                label = ("check:count" if name == "counts" else "check:" + name) if name in CHECKS else name
                self.steps.append([label, c1 - c0 - spare[0], w1 - w0 - spare[1]])
        return call


def kernel_of(leg):
    if leg == "smlib":
        class Bench(cc.Smlib):
            ACCURACIES = (1e-6, 1e-5, 1e-4)
            spare = [0, 0]

            def copy(self, s):
                c0, w0 = time.process_time_ns(), time.perf_counter_ns()
                r = super().copy(s)
                self.spare[0] += time.process_time_ns() - c0
                self.spare[1] += time.perf_counter_ns() - w0
                return r
        return Bench()
    return cc.Occt()


def run_case(leg, case):
    """The build's steps, then each check on a build of its own: what a kernel keeps from one
    check would make the next look quicker than it is."""
    steps = []
    for check in [None] + list(case.get("checks", [])):
        timed = Timed(kernel_of(leg))
        interpreter = cc.Interpreter(timed, {**case, "checks": [check] if check else []})
        records = interpreter.run()
        # (A step the leg cannot do leaves less work timed: such a case does not count.)
        missing = [k for k, v in interpreter.names.items() if v is None]
        missing += [r["measure"] for r in records if r.get("unavailable")]
        if missing:
            raise RuntimeError("not done here: " + ", ".join(missing))
        steps += [s for s in timed.steps if s[0].startswith("check:")] if check else timed.steps
    return steps


def worker(leg, repeats, only, files):
    """The Python kernels' counterpart of the bench_cases binary."""
    if os.environ.get("BENCH_PARALLEL") and leg == "occt":
        from OCP.BOPAlgo import BOPAlgo_Options
        BOPAlgo_Options.SetParallelMode_s(True)
    wanted = set(open(only).read().split()) if only else None
    for path in files:
        stem = os.path.splitext(os.path.basename(path))[0]
        for case in cc.load_cases(path):
            name = f"{stem}::{case['case']}"
            if wanted is not None and name not in wanted:
                continue
            best = None
            try:
                for _ in range(repeats):
                    steps = run_case(leg, case)
                    if best is None:
                        best = steps
                    else:
                        for b, t in zip(best, steps):
                            b[1], b[2] = min(b[1], t[1]), min(b[2], t[2])
                print(json.dumps({"case": name, "steps": best}), flush=True)
            except Exception as e:
                print(json.dumps({"case": name, "error": str(e)[:200]}), flush=True)


def command(leg, args, files, repeats, only):
    pin = ["taskset", "-c", args.pin] if args.pin else []
    if leg == "kentron":
        return pin + [args.kentron] + files
    python = args.smlib_python if leg == "smlib" else args.occt_python
    return pin + [python, os.path.abspath(__file__), "--worker", leg, str(repeats), only or ""] + files


def run_leg(leg, args, files, repeats, only, timeout):
    env = dict(os.environ, REPEATS=str(repeats), OPENBLAS_NUM_THREADS="1")
    if only:
        env["ONLY"] = only
    if args.parallel:
        env["BENCH_PARALLEL"] = "1"
    else:
        env.update(RAYON_NUM_THREADS="1", TBB_NUM_THREADS="1", OMP_NUM_THREADS="1")
    try:
        done = subprocess.run(command(leg, args, files, repeats, only), env=env, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as e:
        return {"timeout": True}, (e.stdout or b"")
    out = {}
    for line in done.stdout.splitlines():
        if line.startswith('{"case"'):
            r = json.loads(line)
            out[r["case"]] = r
    return out, done.returncode


def run(args):
    legs = ("kentron", "occt", "smlib")
    # First each case alone: one that crashes or hangs a kernel is kept out of its rounds.
    names, usable = [], {leg: [] for leg in legs}
    for path in args.files:
        stem = os.path.splitext(os.path.basename(path))[0]
        names += [(path, f"{stem}::{c['case']}") for c in cc.load_cases(path)]
    from concurrent.futures import ThreadPoolExecutor
    tmp = args.out + ".only"
    failures = {leg: {} for leg in legs}

    def probe(job):
        leg, path, name = job
        only = f"{tmp}.{leg}.{abs(hash(name))}"
        with open(only, "w") as f:
            f.write(name + "\n")
        probe_args = argparse.Namespace(**{**vars(args), "pin": None})
        try:
            result, _ = run_leg(leg, probe_args, [path], 1, only, 120)
        finally:
            os.remove(only)
        r = result.get(name) if isinstance(result, dict) else None
        if r and "steps" in r:
            usable[leg].append(name)
        else:
            failures[leg][name] = (r or {}).get("error", "crashed or timed out")
    with ThreadPoolExecutor(10) as pool:
        list(pool.map(probe, [(leg, p, n) for leg in legs for p, n in names]))

    best = {leg: {} for leg in legs}
    for r in range(args.rounds):
        for leg in legs:
            only = f"{tmp}.{leg}"
            with open(only, "w") as f:
                f.write("\n".join(usable[leg]) + "\n")
            result, code = run_leg(leg, args, args.files, 3, only, 3000)
            os.remove(only)
            for name, row in result.items():
                if "steps" not in row:
                    failures[leg].setdefault(name, row.get("error", ""))
                    continue
                have = best[leg].get(name)
                if have is None:
                    best[leg][name] = row["steps"]
                else:
                    for b, t in zip(have, row["steps"]):
                        b[1], b[2] = min(b[1], t[1]), min(b[2], t[2])
            print(f"round {r + 1} {leg}: {len(result)} cases, exit {code}", flush=True)
    with open(args.out, "w") as f:
        json.dump({"best": best, "failures": failures, "parallel": args.parallel, "pin": args.pin}, f)


def report(args):
    data = json.load(open(args.out))
    column = 1 if args.clock == "cpu" else 2
    best, legs = data["best"], ("kentron", "occt", "smlib")
    both = set(best["kentron"]) & set(best["occt"]) & set(best["smlib"])
    def total(leg, case, checks):
        return sum(s[column] for s in best[leg][case] if s[0].startswith("check:") == checks) / 1e6
    def geomean(values):
        values = list(values)
        return math.exp(sum(math.log(v) for v in values) / len(values)) if values else float("nan")
    print(f"clock {args.clock}, pin {data['pin']}, parallel {data['parallel']}; finished: "
          + ", ".join(f"{l} {len(best[l])}" for l in legs) + f"; all three {len(both)}")
    for leg in legs:
        for name, why in sorted(data["failures"][leg].items()):
            print(f"  not timed on {leg}: {name}: {why[:90]}")
    build = {c: {l: total(l, c, False) for l in legs} for c in both}
    print("\nBuilding (ms):        total    geomean/case   median/case")
    for l in legs:
        v = sorted(build[c][l] for c in both)
        print(f"  {l:8} {sum(v):14.1f} {geomean(max(x, 1e-4) for x in v):14.3f} {v[len(v) // 2]:13.3f}")
    floor = 0.02  # Cases of a few microseconds say more about call overheads.
    real = [c for c in both if max(build[c].values()) >= floor]
    print(f"\nRatios over the {len(real)} cases where some kernel took at least {floor} ms (above 1: kentron is faster):")
    for other in ("occt", "smlib"):
        r = sorted(build[c][other] / max(build[c]["kentron"], 1e-4) for c in real)
        wins = sum(x > 1 for x in r)
        print(f"  {other:6}/kentron geomean {geomean(r):7.2f}x  median {r[len(r) // 2]:7.2f}x  kentron faster on {wins}/{len(r)}")
    print("\nBy size of case, by the slowest kernel's time (calls through Python cost microseconds, so"
          "\nthe small ones favour kentron; from 1 ms that does not matter):")
    for low, high in ((0, 0.1), (0.1, 1), (1, 10), (10, 1e9)):
        cs = [c for c in both if low <= max(build[c].values()) < high]
        if cs:
            ro = [build[c]["occt"] / max(build[c]["kentron"], 1e-4) for c in cs]
            rs = [build[c]["smlib"] / max(build[c]["kentron"], 1e-4) for c in cs]
            print(f"  {low:>5}-{high if high < 1e9 else 'up':<5} ms {len(cs):4} cases   occt/kentron {geomean(ro):7.2f}x (kentron faster on {sum(r > 1 for r in ro)})"
                  f"   smlib/kentron {geomean(rs):7.2f}x (kentron faster on {sum(r > 1 for r in rs)})")
    print("\nBy file (ms total, geomean ratio):")
    print(f"  {'file':16}{'cases':>6}{'kentron':>11}{'occt':>11}{'smlib':>11}{'occt/k':>9}{'smlib/k':>9}")
    for stem in sorted({c.split('::')[0] for c in both}):
        cs = [c for c in both if c.startswith(stem + "::")]
        t = {l: sum(build[c][l] for c in cs) for l in legs}
        rs = [c for c in cs if c in real]
        print(f"  {stem:16}{len(cs):6}{t['kentron']:11.2f}{t['occt']:11.2f}{t['smlib']:11.2f}"
              f"{geomean(build[c]['occt'] / max(build[c]['kentron'], 1e-4) for c in rs):9.2f}"
              f"{geomean(build[c]['smlib'] / max(build[c]['kentron'], 1e-4) for c in rs):9.2f}")
    print("\nBy operation (ms per call, geomean over calls in the cases all three finished):")
    ops = {}
    for c in both:
        rows = {l: best[l][c] for l in legs}
        names = {s[0] for l in legs for s in rows[l]}
        for op in names:
            # (Only where all three made the same calls: the operation's times are then paired.)
            calls = {l: [s for s in rows[l] if s[0] == op] for l in legs}
            if len({len(v) for v in calls.values()}) == 1:
                for l in legs:
                    ops.setdefault(op, {}).setdefault(l, []).extend(max(s[column] / 1e6, 1e-4) for s in calls[l])
    print(f"  {'operation':20}{'calls k/o/s':>14}{'kentron':>10}{'occt':>10}{'smlib':>10}")
    for op, by in sorted(ops.items(), key=lambda kv: -sum(sum(v) for v in kv[1].values())):
        # (SMLib's validity and counts read flags where the others check: not the same work.)
        if len(by) == 3 and op not in ("check:valid", "check:count"):
            print(f"  {op:20}{'/'.join(str(len(by[l])) for l in legs):>14}" + "".join(f"{geomean(by[l]):10.3f}" for l in legs))
    print("\nkentron's slowest against the better of the others:")
    order = sorted(real, key=lambda c: build[c]["kentron"] / max(min(build[c]["occt"], build[c]["smlib"]), 1e-4), reverse=True)
    for c in order[:10]:
        print(f"  {c:52}" + "".join(f"{build[c][l]:10.3f}" for l in legs))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5:])
        return
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="action", required=True)
    r = sub.add_parser("run")
    r.add_argument("files", nargs="+")
    r.add_argument("--kentron", required=True)
    r.add_argument("--smlib-python", required=True)
    r.add_argument("--occt-python", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--pin")
    r.add_argument("--parallel", action="store_true")
    r.add_argument("--rounds", type=int, default=5)
    p = sub.add_parser("report")
    p.add_argument("out")
    p.add_argument("--clock", choices=("cpu", "wall"), default="cpu")
    args = parser.parse_args()
    (run if args.action == "run" else report)(args)


if __name__ == "__main__":
    main()
