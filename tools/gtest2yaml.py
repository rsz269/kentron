#!/usr/bin/env python3
"""Moves hand-translated OCCT GTests that only build shapes and check measurable facts from
Rust (crates/occt-tests/tests/gtest_*) to YAML cases (docs/test-cases.md).

Each test body is parsed as a sequence of the fixture helpers' calls (create_box,
perform_union, validate_result, ...). A test whose every statement maps onto the case format
becomes a YAML case in crates/occt-tests/cases/occt/gtest/<module>/<file>.yaml, with
ported_from naming the OCCT C++ file and test; it is removed from the Rust file. Others stay.

Usage:
    tools/gtest2yaml.py tkbo [--occt ~/OCCT] [--dry-run]
"""

import argparse
import json
import math
import re
import subprocess
from pathlib import Path

import drawyaml
from drawyaml import Expr, Ref, flow, key

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "crates" / "occt-tests" / "tests"
CASES = ROOT / "crates" / "occt-tests" / "cases" / "occt" / "gtest"

SOURCES = {"tkbo": "src/ModelingAlgorithms/TKBO/GTests", "tkprim": "src/ModelingAlgorithms/TKPrim/GTests"}


class Unsupported(Exception):
    pass


# --- Rust expressions ----------------------------------------------------------------------

TOKEN = re.compile(r"\s*(?:(\d+\.\d*(?:e[-+]?\d+)?|\d+(?:e[-+]?\d+)?)(?:_?f64|_?usize)?|([A-Za-z_][A-Za-z0-9_]*(?:::[A-Za-z_][A-Za-z0-9_]*)*!?)|(&mut|\.\.|[-+*/(),\[\]&.;|<>=]))")


def tokens(text):
    out, at = [], 0
    while at < len(text):
        if text[at:].strip() == "":
            break
        m = TOKEN.match(text, at)
        if not m or m.end() == at:
            raise Unsupported(f"token at {text[at:at + 20]!r}")
        out.append(m.group(1) or m.group(2) or m.group(3))
        at = m.end()
    return out


class Parser:
    def __init__(self, toks):
        self.t, self.i = toks, 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None

    def take(self, expected=None):
        tok = self.peek()
        if tok is None or (expected is not None and tok != expected):
            raise Unsupported(f"expected {expected}, got {tok}")
        self.i += 1
        return tok

    def expr(self):
        v = self.term()
        while self.peek() in ("+", "-"):
            op = self.take()
            v = ("bin", op, v, self.term())
        return v

    def term(self):
        v = self.unary()
        while self.peek() in ("*", "/"):
            op = self.take()
            v = ("bin", op, v, self.unary())
        return v

    def unary(self):
        if self.peek() == "-":
            self.take()
            return ("neg", self.unary())
        if self.peek() in ("&", "&mut"):
            self.take()
            return self.unary()
        return self.postfix()

    def postfix(self):
        v = self.atom()
        while self.peek() == ".":
            self.take()
            name = self.take()
            args = self.args() if self.peek() == "(" else None
            v = ("method", v, name, args)
        return v

    def args(self):
        self.take("(")
        out = []
        while self.peek() != ")":
            out.append(self.expr())
            if self.peek() == ",":
                self.take()
        self.take(")")
        return out

    def atom(self):
        tok = self.take()
        if tok == "(":
            items = [self.expr()]
            while self.peek() == ",":
                self.take()
                if self.peek() == ")":
                    break
                items.append(self.expr())
            self.take(")")
            return items[0] if len(items) == 1 else ("tuple", items)
        if tok == "[":
            items = []
            while self.peek() != "]":
                items.append(self.expr())
                if self.peek() == ",":
                    self.take()
            self.take("]")
            return ("list", items)
        if tok == "vec!":
            return self.atom()
        if re.match(r"\d", tok):
            return ("num", float(tok))
        if self.peek() == "(":
            return ("call", tok, self.args())
        return ("name", tok)


# --- values --------------------------------------------------------------------------------

class Conv:
    def __init__(self):
        self.env = {}  # Rust local -> value
        self.steps = []
        self.tmp = 0
        self.boxes = {}  # shape name -> (corner, size)
        self.unique_of = None

    def fresh(self):
        self.tmp += 1
        return f"t{self.tmp}"

    def bind(self, name, op, box=None):
        self.steps.append(("bind", name, op))
        if box:
            self.boxes[name] = box
        return ("shape", name)

    def shape_op(self, op, name=None, box=None):
        return self.bind(name or self.fresh(), op, box)

    def num(self, e):
        v = self.value(e)
        if v[0] != "num":
            raise Unsupported(f"not a number: {e}")
        return v[1]

    def shape(self, e):
        v = self.value(e)
        if v[0] != "shape":
            raise Unsupported(f"not a shape: {e}")
        return Ref(v[1])

    def vec(self, e):
        v = self.value(e)
        if v[0] not in ("pt", "vec", "dir"):
            raise Unsupported(f"not a vector: {e}")
        return v[1]

    def value(self, e):
        kind = e[0]
        if kind == "num":
            return e
        if kind == "neg":
            v = self.value(e[1])
            if v[0] == "num":
                return ("num", -v[1])
            if v[0] in ("dir", "vec", "pt"):
                return (v[0], [-x for x in v[1]])
            raise Unsupported("negation")
        if kind == "bin":
            a, b = self.value(e[2]), self.value(e[3])
            if a[0] == b[0] == "num":
                return ("num", {"+": a[1] + b[1], "-": a[1] - b[1], "*": a[1] * b[1], "/": a[1] / b[1]}[e[1]])
            if a[0] == "pt" and b[0] == "vec" and e[1] in "+-":
                s = 1 if e[1] == "+" else -1
                return ("pt", [x + s * y for x, y in zip(a[1], b[1])])
            raise Unsupported("arithmetic")
        if kind == "name":
            n = e[1]
            consts = {"PI": math.pi, "SQRT_2": math.sqrt(2), "std::f64::consts::PI": math.pi,
                      "kentron::tolerance::LINEAR": 1e-7, "DEFAULT_TOLERANCE": 1e-6}
            if n in consts:
                return ("num", consts[n])
            dirs = {"Direction::X": [1, 0, 0], "Direction::Y": [0, 1, 0], "Direction::Z": [0, 0, 1]}
            if n in dirs:
                return ("dir", dirs[n])
            axes = {"Axis::X": [1, 0, 0], "Axis::Y": [0, 1, 0], "Axis::Z": [0, 0, 1]}
            if n in axes:
                return ("axis", [0, 0, 0], axes[n])
            if n == "Point::ORIGIN":
                return ("pt", [0, 0, 0])
            if n in ("true", "false"):
                return ("bool", n == "true")
            if n.startswith("NodeKind::"):
                return ("kind", n[len("NodeKind::"):])
            if n in self.env:
                return self.env[n]
            raise Unsupported(f"name {n}")
        if kind == "list":
            return ("list", [self.value(x) for x in e[1]])
        if kind == "tuple":
            return ("tuple", [self.value(x) for x in e[1]])
        if kind == "method":
            return self.method(e)
        if kind == "call":
            return self.call(e[1], e[2])
        raise Unsupported(f"expression {kind}")

    def method(self, e):
        _, obj, name, args = e
        if name == "moved" and args:
            s = self.shape(obj)
            t = self.value(args[0])
            if t[0] == "translation":
                return self.shape_op({"translate": {"shape": s, "by": t[1]}})
            if t[0] == "rotation":
                return self.shape_op({"rotate": {"shape": s, "axis": {"origin": t[1], "direction": t[2]},
                                                 "angle": t[3] * 180 / math.pi}})
            raise Unsupported("moved")
        if name == "sqrt" and args == []:
            return ("num", math.sqrt(self.num(obj)))
        raise Unsupported(f"method {name}")

    def call(self, f, args):
        a = args
        if f == "pnt":
            return ("pt", [self.num(x) for x in a])
        if f in ("vec3", "Vector::new"):
            return ("vec", [self.num(x) for x in a])
        if f in ("dir", "Direction::new"):
            return ("dir", [self.num(x) for x in a])
        if f == "Point::new":
            return ("pt", [self.num(x) for x in a])
        if f == "Axis::new":
            return ("axis", self.vec(a[0]), self.vec(a[1]))
        if f == "Frame::new":
            return ("frame", {"origin": self.vec(a[0]), "z": self.vec(a[1]), "x": self.vec(a[2])})
        if f == "Plane::new":
            return ("frame", {"origin": self.vec(a[0]), "z": self.vec(a[1])})
        if f == "Plane::from_frame":
            v = self.value(a[0])
            if v[0] != "frame":
                raise Unsupported("plane")
            return v
        if f == "Transform::translation":
            return ("translation", self.vec(a[0]))
        if f == "Transform::rotation":
            ax = self.value(a[0])
            return ("rotation", ax[1], ax[2], self.num(a[1]))
        if f == "op":
            cmd = a[0]
            if cmd[0] != "name" or not cmd[1].startswith("Cmd::"):
                raise Unsupported("op")
            params = self.value(a[1])
            return ("op", [cmd[1][len("Cmd::"):]] + [p[1] for p in params[1]])
        # shapes
        if f == "create_unit_sphere":
            return self.shape_op({"sphere": {"radius": 1}})
        if f == "create_unit_box":
            return self.shape_op({"box": {"size": [1, 1, 1]}}, box=([0, 0, 0], [1, 1, 1]))
        if f == "create_box":
            corner, size = self.vec(a[1]), [self.num(x) for x in a[2:5]]
            return self.shape_op({"box": {"corner": corner, "size": size}}, box=(corner, size))
        if f == "create_identical_boxes":
            return ("tuple", [self.call("create_unit_box", a), self.call("create_unit_box", a)])
        if f == "create_sphere":
            return self.shape_op({"sphere": {"radius": self.num(a[2]), "frame": {"origin": self.vec(a[1]), "z": [0, 0, 1], "x": [1, 0, 0]}}})
        if f == "create_cylinder":
            return self.shape_op({"cylinder": {"radius": self.num(a[1]), "height": self.num(a[2])}})
        if f == "create_cone":
            return self.shape_op({"cone": {"r1": self.num(a[1]), "r2": self.num(a[2]), "height": self.num(a[3])}})
        if f == "create_cylinder_on_plane":
            frame = self.value(a[1])
            return self.shape_op({"cylinder": {"radius": self.num(a[2]), "height": self.num(a[3]), "frame": frame[1]}})
        if f == "convert_to_nurbs":
            return self.shape_op({"to_nurbs": {"shape": self.shape(a[1])}})
        if f in ("rotate_shape",):
            ax = self.value(a[1])
            return self.shape_op({"rotate": {"shape": self.shape(a[0]), "axis": {"origin": ax[1], "direction": ax[2]},
                                             "angle": self.num(a[2]) * 180 / math.pi}})
        if f in ("rotate_x", "rotate_y", "rotate_z"):
            d = {"rotate_x": [1, 0, 0], "rotate_y": [0, 1, 0], "rotate_z": [0, 0, 1]}[f]
            return self.shape_op({"rotate": {"shape": self.shape(a[0]), "axis": {"origin": [0, 0, 0], "direction": d},
                                             "angle": self.num(a[1])}})
        if f == "rotate_standard":
            s = self.call("rotate_z", [a[0], ("num", -90.0)])
            self.env["__std"] = s
            return self.call("rotate_y", [("name", "__std"), ("num", -45.0)])
        if f == "rotate_y90":
            return self.shape_op({"rotate": {"shape": self.shape(a[0]), "axis": {"origin": [0, 0, 1], "direction": [0, 1, 0]},
                                             "angle": 90}})
        if f == "translate_shape":
            return self.shape_op({"translate": {"shape": self.shape(a[0]), "by": self.vec(a[1])}})
        if f == "translate":
            return self.shape_op({"translate": {"shape": self.shape(a[0]), "by": [self.num(x) for x in a[1:4]]}})
        if f == "create_vertex":
            return self.shape_op({"vertex": {"at": self.vec(a[1])}})
        if f == "create_edge":
            return self.shape_op({"edge": {"from": self.vec(a[1]), "to": self.vec(a[2])}})
        if f == "create_polygon_wire":
            pts = [p[1] for p in self.value(a[1])[1]]
            close = self.value(a[2])
            if close[0] != "bool":
                raise Unsupported("close")
            if close[1]:
                pts = pts + [pts[0]]
            return self.shape_op({"polyline": {"points": pts}})
        if f == "create_polygon_face":
            pts = [self.vec(x) for x in a[1:5]]
            w = self.shape_op({"polyline": {"points": pts + [pts[0]]}})
            return self.shape_op({"face": {"wire": Ref(w[1])}})
        if f == "create_rectangular_profile":
            c, x, y = self.vec(a[1]), self.num(a[2]), self.num(a[3])
            pts = [c, [c[0] + x, c[1], c[2]], [c[0] + x, c[1] + y, c[2]], [c[0], c[1] + y, c[2]], c]
            return self.shape_op({"polyline": {"points": pts}})
        if f == "create_rectangular_prism":
            w = self.call("create_rectangular_profile", a[:4])
            face = self.shape_op({"face": {"wire": Ref(w[1])}})
            return self.shape_op({"extrude": {"profile": Ref(face[1]), "by": [0, 0, self.num(a[4])]}})
        if f == "create_face_from_wire":
            return self.shape_op({"face": {"wire": self.shape(a[1])}})
        if f == "create_prism":
            return self.shape_op({"extrude": {"profile": self.shape(a[1]), "by": self.vec(a[2])}})
        if f == "create_revolution":
            ax = self.value(a[2])
            return self.shape_op({"revolve": {"profile": self.shape(a[1]), "axis": {"origin": ax[1], "direction": ax[2]},
                                              "angle": self.num(a[3]) * 180 / math.pi}})
        if f in ("create_profile", "create_profile_from_operations"):
            if f == "create_profile":
                frame, ops = self.value(a[1]), self.value(a[2])
                frame = frame[1]
            else:
                frame, ops = None, self.value(a[1])
            op = drawyaml.readable_profile([o[1] for o in ops[1]])
            if frame is not None:
                op = {"frame": frame, **op}
            return self.shape_op({"profile": op})
        if f == "create_blend":
            s = self.shape(a[1])
            if self.num(a[2]) != 1 or s not in self.boxes:
                raise Unsupported("create_blend on other than a box's first edge")
            corner, size = self.boxes[s]
            # The first edge OCCT's explorer meets on a box: the vertical one at its corner.
            e = self.shape_op({"select": {"of": s, "edge": [corner[0], corner[1], corner[2] + size[2] / 2]}})
            return self.shape_op({"fillet": {"shape": s, "edges": [{"edge": Ref(e[1]), "radius": self.num(a[3])}]}})
        if f in ("perform_union", "perform_difference", "perform_intersection"):
            op = f[len("perform_"):]
            return self.shape_op({op: {"object": self.shape(a[1]), "tool": self.shape(a[2])}})
        raise Unsupported(f"call {f}")

    # --- statements ------------------------------------------------------------------------

    def statement(self, text):
        text = text.strip()
        if not text or text in ("let mut model = Model::new()",):
            return
        m = re.match(r"let unique = \|k(?:: NodeKind)?\| model\.unique_subshapes\((\w+), k\)\.len\(\)$", text)
        if m:
            self.unique_of = m.group(1)
            return
        m = re.match(r"let (?:mut )?\((\w+), (\w+)\) = (.*)$", text, re.S)
        if m:
            v = self.value(Parser(tokens(m.group(3))).expr())
            if v[0] != "tuple":
                raise Unsupported("tuple")
            for n, x in zip(m.group(1, 2), v[1]):
                self.assign(n, x)
            return
        m = re.match(r"(?:let (?:mut )?)?(\w+) = (.*)$", text, re.S)
        if m and not text.startswith(("assert", "validate_result")):
            self.assign(m.group(1), self.value(Parser(tokens(m.group(2))).expr()))
            return
        self.assertion(text)

    def assign(self, name, v):
        if v[0] == "shape":
            # Bind under the Rust name (rename the last step if it made this value).
            last = self.steps[-1] if self.steps else None
            if last and last[0] == "bind" and last[1] == v[1] and v[1].startswith("t"):
                self.steps[-1] = ("bind", name, last[2])
                if v[1] in self.boxes:
                    self.boxes[name] = self.boxes.pop(v[1])
                v = ("shape", name)
            else:
                self.steps.append(("bind", name, Ref(v[1])))
                if v[1] in self.boxes:
                    self.boxes[name] = self.boxes[v[1]]
                v = ("shape", name)
        self.env[name] = v

    def check(self, c):
        self.steps.append(("check", c))

    def assertion(self, text):
        m = re.match(r"validate_result\(&model, (\w+), (.*)\)$", text, re.S)
        if m:
            name = self.shape(("name", m.group(1)))
            p = Parser(tokens(m.group(2)))
            area = self.num(p.expr())
            p.take(",")
            volume = self.num(p.expr())
            p.take(",")
            empty = self.value(p.expr())
            if empty[0] != "bool":
                raise Unsupported("validate_result")
            if empty[1]:
                self.check({"area": name, "le": 1e-6})
                return
            if area >= 0:
                self.check({"area": name, "eq": area, "abs": 1e-4 * max(area, 1.0)})
            if volume >= 0:
                self.check({"volume": name, "eq": volume, "abs": 1e-6})
            return
        m = re.match(r"assert_near!\(get_(surface_area|volume)\(&model, (\w+)\), (.*)\)$", text, re.S)
        if m:
            name = self.shape(("name", m.group(2)))
            p = Parser(tokens(m.group(3)))
            v = self.num(p.expr())
            p.take(",")
            tol = self.num(p.expr())
            self.check({"area" if m.group(1) == "surface_area" else "volume": name, "eq": v, "abs": tol})
            return
        m = re.match(r"assert_eq!\(model\.explore\((\w+), NodeKind::(\w+)\)\.count\(\), (\d+)\)$", text)
        if m:
            self.count(m.group(1), m.group(2), int(m.group(3)))
            return
        m = re.match(r"assert_eq!\(unique\(NodeKind::(\w+)\), (\d+)\)$", text)
        if m and self.unique_of:
            self.count(self.unique_of, m.group(1), int(m.group(2)))
            return
        m = re.match(r"assert!\(validity::check\(&model, (\w+)\)\.is_valid\(\)\)$", text)
        if m:
            self.check({"valid": self.shape(("name", m.group(1)))})
            return
        m = re.match(r"assert_eq!\((\w+)\.kind\(\), NodeKind::(\w+)\)$", text)
        if m:
            self.check({"kind": self.shape(("name", m.group(1))), "is": m.group(2).lower()})
            return
        raise Unsupported(f"statement {text[:50]!r}")

    def count(self, name, kind, n):
        plural = {"Vertex": "vertices", "Edge": "edges", "Wire": "wires", "Face": "faces", "Shell": "shells",
                  "Solid": "solids", "CompSolid": "compsolids", "Compound": "compounds"}[kind]
        self.check({"count": self.shape(("name", name)), "eq": {plural: n}})


# --- files ---------------------------------------------------------------------------------

def split_statements(body):
    out, depth, cur = [], 0, ""
    for ch in body:
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == ";" and depth == 0:
            out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur.strip():
        out.append(cur)
    return out


def tests_of(text):
    """(start, end, name, body, comments) of each top-level #[test] fn."""
    out = []
    for m in re.finditer(r"((?:///[^\n]*\n)*)#\[test\]\n(?:#\[[^\n]*\]\n)*fn (\w+)\(\) \{\n", text):
        start, at, depth = m.start(), m.end(), 1
        while depth:
            depth += {"{": 1, "}": -1}.get(text[at], 0)
            at += 1
        out.append((start, at, m.group(2), text[m.end():at - 1], m.group(1)))
    return out


def snake(name):
    s = re.sub(r"(?<=[a-z0-9])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])", "_", name)
    return s.lower()


def occt_names(occt, module):
    """Rust-style name -> (C++ file, Fixture.Test)."""
    out = {}
    for f in sorted((occt / SOURCES[module]).glob("*.cxx")):
        for m in re.finditer(r"^TEST(?:_F)?\((\w+),\s*(\w+)\)", f.read_text(), re.M):
            out.setdefault(snake(m.group(2)), (f.relative_to(occt).as_posix(), f"{m.group(1)}.{m.group(2)}"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("module")
    ap.add_argument("--occt", default=str(Path.home() / "OCCT"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    occt = Path(args.occt).expanduser()
    git_hash = subprocess.run(["git", "-C", str(occt), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    names = occt_names(occt, args.module)
    moved = kept = 0
    for rs in sorted((TESTS / f"gtest_{args.module}").glob("*.rs")):
        if rs.name == "main.rs":
            continue
        text = rs.read_text()
        docs, remove = [], []
        for start, end, name, body, comments in tests_of(text):
            origin = names.get(name)
            if origin is None:
                kept += 1
                continue
            conv = Conv()
            notes = [c.strip()[3:].strip() for c in comments.splitlines() if c.strip()]
            body_lines = []
            for line in body.splitlines():
                if line.strip().startswith("//"):
                    notes.append(line.strip()[2:].strip())
                else:
                    body_lines.append(line)
            try:
                for stmt in split_statements("\n".join(body_lines)):
                    conv.statement(" ".join(stmt.split()))
            except (Unsupported, drawyaml.Unsupported, IndexError, KeyError, TypeError) as e:
                kept += 1
                continue
            docs.append(document(name, origin, git_hash, conv, notes))
            remove.append((start, end))
            moved += 1
        if not docs:
            continue
        out = CASES / args.module / rs.with_suffix(".yaml").name
        if not args.dry_run:
            header = f"# OCCT GTests moved from crates/occt-tests/tests/gtest_{args.module}/{rs.name} by tools/gtest2yaml.py.\n"
            existing = out.read_text() if out.exists() else header
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(existing + ("---\n" if existing != header else "") + "---\n".join(docs))
            for start, end in reversed(remove):
                text = text[:start] + text[end:].lstrip("\n")
            rs.write_text(text)
    print(f"moved {moved}, kept {kept}")


def document(name, origin, git_hash, conv, notes):
    last_bind = max((i for i, s in enumerate(conv.steps) if s[0] == "bind"), default=-1)
    build, checks = conv.steps[:last_bind + 1], [s[1] for s in conv.steps[last_bind + 1:]]
    out = [f"case: {key(name)}", "ported_from:", "  project: OCCT", f"  repo: {drawyaml.OCCT_REPO}",
           f"  git_hash: {git_hash}", f"  file: {origin[0]}", f"  test: {origin[1]}"]
    if notes:
        out.append("notes:")
        out += [f"  - {json.dumps(n)}" for n in notes]
    out.append("build:")
    for s in build:
        out.append(f"  - {key(s[1])}: {flow(s[2])}" if s[0] == "bind" else f"  - check: {flow(s[1])}")
    if checks:
        out.append("checks:")
        out += [f"  - {flow(c)}" for c in checks]
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    main()
