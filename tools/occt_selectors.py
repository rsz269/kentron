#!/usr/bin/env python3
"""Pins down the sub-shapes that OCCT DRAW scripts name by enumeration order.

A script like `explode s E; blend result s 10 s_3` picks "the third edge OCCT's explorer
visits". Kentron enumerates differently, so the ported test would fillet another edge. This
tool runs each such script in OCCT's DRAW (a local build), records the centroid of every
sub-shape an `explode`/`nexplode` creates (of length for edges, of area for faces, the point
for vertices), and writes those the script uses afterwards to
crates/occt-tests/draw_selectors.json. tools/draw2rs.py then rebinds the names to the
sub-shapes at those positions (kentron's `nearest` command).

Usage:
    tools/occt_selectors.py [--occt ~/OCCT] [--draw DIR] [--tcl-lib DIR] [--groups blend,...]

`--draw` is DRAW's install directory (default `$OCCT_DRAW`, else `~/occt-build/install`);
`--tcl-lib` the Tcl libraries DRAW was built with, if not the system's (default
`$OCCT_TCL_LIB`, else `tcl/lib` beside the install directory if there is one).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "crates" / "occt-tests" / "draw_selectors.json"
# OCCT's test dataset (an external download, see crates/occt-tests/README.md).
DATA = os.environ.get("OCCT_TEST_DATA", str(Path.home() / "occt-data"))

EXPLODE = re.compile(r"^\s*(n?explode)\s+(\S+)\s+([EeFfVvWw]|[Ss][Oo]|[Ss][Hh])\b")

WRAPPER = r"""
set ::kentron_k 0
proc kentron_report {k n t} {
  switch -nocase -- $t {
    e {
      # Sampled along the curve: lprops is imprecise on some B-spline edges.
      if {[catch {uplevel #0 [list mkcurve kentron_c $n]}]} { return }
      uplevel #0 {bounds kentron_c kentron_u1 kentron_u2}
      set u1 [uplevel #0 {dval kentron_u1}]; set u2 [uplevel #0 {dval kentron_u2}]
      set m 4000; set len 0.; set sx 0.; set sy 0.; set sz 0.
      uplevel #0 [list cvalue kentron_c $u1 kentron_x kentron_y kentron_z]
      set px [uplevel #0 {dval kentron_x}]; set py [uplevel #0 {dval kentron_y}]; set pz [uplevel #0 {dval kentron_z}]
      for {set i 1} {$i <= $m} {incr i} {
        uplevel #0 [list cvalue kentron_c [expr {$u1 + ($u2 - $u1) * $i / double($m)}] kentron_x kentron_y kentron_z]
        set qx [uplevel #0 {dval kentron_x}]; set qy [uplevel #0 {dval kentron_y}]; set qz [uplevel #0 {dval kentron_z}]
        set d [expr {sqrt(($qx-$px)**2 + ($qy-$py)**2 + ($qz-$pz)**2)}]
        set len [expr {$len + $d}]
        set sx [expr {$sx + $d * ($px + $qx) / 2.}]; set sy [expr {$sy + $d * ($py + $qy) / 2.}]; set sz [expr {$sz + $d * ($pz + $qz) / 2.}]
        set px $qx; set py $qy; set pz $qz
      }
      if {$len > 0} { puts [format "KENTRON_SEL %d %s E %.12g %.12g %.12g %.12g" $k $n [expr {$sx/$len}] [expr {$sy/$len}] [expr {$sz/$len}] $len] }
      return
    }
    f { set p [uplevel #0 [list sprops $n]] }
    sh { set p [uplevel #0 [list sprops $n]] }
    so { set p [uplevel #0 [list vprops $n]] }
    w { set p [uplevel #0 [list lprops $n]] }
    v {
      set b [regexp -all -inline {[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?} [uplevel #0 [list bounding $n]]]
      puts "KENTRON_SEL $k $n V 0 [expr {([lindex $b 0] + [lindex $b 3]) / 2.}] [expr {([lindex $b 1] + [lindex $b 4]) / 2.}] [expr {([lindex $b 2] + [lindex $b 5]) / 2.}]"
      return
    }
  }
  if {[regexp {Mass : +([-+0-9.eE]+)} $p -> m] && [regexp {X = +([-+0-9.eE]+)} $p -> x] && [regexp {Y = +([-+0-9.eE]+)} $p -> y] && [regexp {Z = +([-+0-9.eE]+)} $p -> z]} {
    puts "KENTRON_SEL $k $n [string toupper $t] $x $y $z $m"
  }
}
foreach kentron_cmd {explode nexplode} {
  rename $kentron_cmd kentron_orig_$kentron_cmd
  proc $kentron_cmd {args} [string map [list CMD $kentron_cmd] {
    # At global level: DRAW creates its variables in the calling scope.
    set r [uplevel #0 [list kentron_orig_CMD {*}$args]]
    set t [lindex $args 1]
    if {[regexp -nocase {^([efvw]|so|sh)$} $t]} {
      foreach n $r { kentron_report $::kentron_k $n $t }
    }
    incr ::kentron_k
    return $r
  }]
}
"""


def explodes(text: str):
    """The script's explode/nexplode commands in order: (shape, type, line index)."""
    out = []
    for i, line in enumerate(text.splitlines()):
        for part in line.split(";"):
            m = EXPLODE.match(part)
            # tests/offset's OFFSETSHAPE explodes `s` into faces when it is given some.
            offset = re.match(r"^\s*OFFSETSHAPE\s+\S+\s+\{\s*\S", part)
            if m:
                out.append((m.group(2), m.group(3).upper(), i))
            elif offset:
                out.append(("s", "F", i))
            elif re.match(r"^\s*n?explode\b", part):
                out.append((None, None, i))
    return out


def used_names(text: str, found: list) -> list:
    """Per explode: the names it creates that the script uses before re-exploding the shape."""
    lines = text.splitlines()
    result = []
    for k, (shape, kind, line) in enumerate(found):
        if shape is None:
            result.append(set())
            continue
        end = next((l for (s, _, l) in found[k + 1:] if s == shape), len(lines))
        # OFFSETSHAPE uses the faces it explodes on its own line.
        first = line if re.match(r"^\s*OFFSETSHAPE\b", lines[line]) else line + 1
        rest = "\n".join(lines[first:end])
        result.append(set(re.findall(rf"\b{re.escape(shape)}_\d+\b", rest)))
    return result


def run_draw(draw: Path, tcl_lib: Path | None, script: str) -> str:
    with tempfile.NamedTemporaryFile("w", suffix=".tcl", delete=False) as f:
        f.write(script)
        path = f.name
    libraries = f'export LD_LIBRARY_PATH="{tcl_lib}:$LD_LIBRARY_PATH"; ' if tcl_lib else ""
    command = (f'{libraries}source "{draw}/bin/env.sh" >/dev/null 2>&1; '
               f'export CSF_TestDataPath="{DATA}"; "{draw}/bin/DRAWEXE" -b -f "{path}"')
    try:
        r = subprocess.run(["bash", "-c", command], capture_output=True, text=True, timeout=180)
        return r.stdout + r.stderr
    except subprocess.TimeoutExpired as e:
        out = e.stdout or ""
        return out.decode() if isinstance(out, bytes) else out
    finally:
        os.unlink(path)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--occt", default=str(Path.home() / "OCCT"))
    ap.add_argument("--draw", default=os.environ.get("OCCT_DRAW", str(Path.home() / "occt-build" / "install")))
    ap.add_argument("--tcl-lib", default=os.environ.get("OCCT_TCL_LIB"))
    ap.add_argument("--groups", default="blend")
    args = ap.parse_args()
    occt, draw = Path(args.occt).expanduser(), Path(args.draw).expanduser()
    beside = draw.parent / "tcl" / "lib"
    tcl_lib = Path(args.tcl_lib).expanduser() if args.tcl_lib else beside if beside.is_dir() else None
    data = json.loads(OUT.read_text()) if OUT.exists() else {}
    data["_comment"] = ("Sub-shapes that OCCT scripts name by enumeration order, by position "
                        "(written by tools/occt_selectors.py from OCCT's DRAW; see tools/draw2rs.py).")
    for group in args.groups.split(","):
        gdir = occt / "tests" / group
        begin = [gdir / "begin"]
        for grid in sorted(p for p in gdir.iterdir() if p.is_dir()):
            for case in sorted(p for p in grid.iterdir() if p.is_file() and p.name not in ("begin", "end")):
                text = case.read_text(errors="replace")
                found = explodes(text)
                names = used_names(text, found)
                if not any(names):
                    continue
                prelude = "pload ALL\n" + "".join(b.read_text() for b in begin + [grid / "begin"] if b.exists())
                output = run_draw(draw, tcl_lib, prelude + WRAPPER + text + "\nexit\n")
                positions: dict = {}
                for line in output.splitlines():
                    m = re.match(r"KENTRON_SEL (\d+) (\S+) (E|F|V|W|SO|SH) (.*)$", line.strip())
                    if m:
                        values = [float(x) for x in m.group(4).split()]
                        if m.group(3) == "V":
                            values = values[1:] + [0.0]
                        positions[(int(m.group(1)), m.group(2))] = [m.group(3)] + values
                entries = []
                for k, wanted in enumerate(names):
                    picked = {n: positions[(k, n)] for n in sorted(wanted) if (k, n) in positions}
                    if picked:
                        entries.append({"explode": k, "names": picked})
                rel = case.relative_to(occt).as_posix()
                if entries:
                    data[rel] = entries
                    print(f"{rel}: {sum(len(e['names']) for e in entries)} names")
                else:
                    print(f"{rel}: nothing recorded")
    OUT.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
