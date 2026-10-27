"""DRAW scripts as YAML test cases (docs/test-cases.md).

tools/draw2rs.py translates a script into DRAW-session calls (`d.box(&["b", "1", "2", "3"])?;`),
with kentron's deviations and position selectors applied. This module turns such a call
sequence into a YAML case when every command maps onto the language; otherwise it raises
`Unsupported` and the script stays a generated Rust test.
"""

from __future__ import annotations

import json
import re

OCCT_REPO = "https://github.com/Open-Cascade-SAS/OCCT"


class Unsupported(Exception):
    pass


CALL = re.compile(r'^d\.(?:r#)?(\w+)\(&\[(.*)\]\)\?;$')
REQUIRE = re.compile(r'^d\.require_var\(("(?:[^"\\]|\\.)*")\)\?;$')
DATA_FILE = re.compile(r'data_file\(("(?:[^"\\]|\\.)*")\)\.as_str\(\)')


def parse_calls(body: list) -> list:
    """(command, [args]) per statement; the arguments must all be constant strings."""
    calls = []
    for stmt in body:
        r = REQUIRE.match(stmt.strip())
        if r:
            calls.append(("require_var", [json.loads(r.group(1))]))
            continue
        m = CALL.match(stmt.strip())
        if not m:
            raise Unsupported(f"statement: {stmt[:40]}")
        inner = DATA_FILE.sub(lambda d: json.dumps("@data:" + json.loads(d.group(1))), m.group(2))
        try:
            args = json.loads("[" + inner + "]") if inner.strip() else []
        except json.JSONDecodeError:
            raise Unsupported(f"arguments: {stmt[:40]}")
        if not all(isinstance(a, str) for a in args):
            raise Unsupported("arguments")
        calls.append((m.group(1), args))
    return calls


NUMBER = re.compile(r"^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")
NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def num(a: str, params: dict):
    """A number argument: a literal, or an expression over the parameters (kept as text)."""
    if NUMBER.match(a):
        v = float(a)
        return int(v) if v.is_integer() and "." not in a and "e" not in a.lower() else v
    names = re.findall(r"[A-Za-z_][A-Za-z0-9_]*", a)
    allowed = {"pi", "PI", "sqrt", "sin", "cos", "tan", "asin", "acos", "atan", "abs", "exp", "log", "atan2",
               "pow", "min", "max"}
    if not names or any(n not in params and n not in allowed for n in names):
        raise Unsupported(f"number {a!r}")
    if not re.fullmatch(r"[A-Za-z0-9_.+\-*/^(), ]+", a):
        raise Unsupported(f"number {a!r}")
    return Expr(a)


class Expr(str):
    """An expression string (always quoted)."""


def substitute(text: str, params: dict, depth: int = 0) -> str:
    """An expression with the parameters' values written in (parenthesized)."""
    def value(m):
        name = m.group(0)
        if name not in params or depth > 20:
            return name
        v = params[name]
        return f"({substitute(v, params, depth + 1)})" if isinstance(v, Expr) else f"({v})"
    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", value, text)


def inline(x, params: dict):
    """`x` (steps, operations, values) with expressions over parameters inlined: evaluated
    when they are arithmetic on numbers, else kept as expressions of numbers."""
    if isinstance(x, Expr):
        text = substitute(x, params)
        if re.fullmatch(r"[0-9.eE+\-*/^() ]+", text):
            v = eval(text.replace("^", "**"), {"__builtins__": {}})
            v = float(f"{v:.12g}")
            return int(v) if v.is_integer() and abs(v) < 1e15 else v
        # Drop the parentheses around plain numbers.
        text = re.sub(r"(?<![A-Za-z0-9_])\((-?[0-9.]+(?:[eE][-+]?[0-9]+)?)\)", r"\1", text)
        # Parentheses around the whole expression.
        while text.startswith("(") and text.endswith(")"):
            depth = 0
            for i, c in enumerate(text):
                depth += (c == "(") - (c == ")")
                if depth == 0 and i < len(text) - 1:
                    break
            else:
                text = text[1:-1]
                continue
            break
        return Expr(text)
    if isinstance(x, Ref):
        return x
    if isinstance(x, dict):
        return {k: inline(v, params) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return type(x)(inline(v, params) for v in x)
    return x


class Ref(str):
    """A name of a bound value."""


def flow(v) -> str:
    """YAML flow syntax."""
    if isinstance(v, dict):
        return "{" + ", ".join(f"{key(k)}: {flow(x)}" for k, x in v.items()) + "}"
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(flow(x) for x in v) + "]"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, Ref) and NAME.match(v) and v.lower() not in RESERVED:
        return str(v)
    return json.dumps(str(v))


RESERVED = {"true", "false", "null", "yes", "no", "on", "off", "y", "n", "~"}


def key(k: str) -> str:
    return k if NAME.match(k) and k.lower() not in RESERVED else json.dumps(k)


# DRAW's profile codes: the number of values each takes.
PROFILE_ARGS = {"F": 2, "O": 3, "P": 6, "S": 1, "X": 1, "XX": 1, "Y": 1, "YY": 1, "L": 1, "T": 2, "TT": 2,
                "R": 1, "RR": 1, "D": 2, "C": 2, "IX": 1, "IY": 1, "W": 0, "WW": 0}
# DRAW's profile codes as the YAML steps (docs/test-cases.md).
PROFILE_STEPS = {"F": "start", "O": "origin", "S": "on", "X": "along_x", "XX": "to_x", "Y": "along_y",
                 "YY": "to_y", "L": "forward", "T": "by", "TT": "to", "R": "turn", "RR": "heading",
                 "D": "direction", "IX": "forward_to_x", "IY": "forward_to_y"}


def profile_code(word: str) -> str:
    """DRAW reads a profile code's first letter, and the second for XX, YY, TT, RR, IX, IY, WW."""
    w = word.upper()
    if len(w) > 1 and w[0] in "XYTRW" and w[1] == w[0]:
        return w[:2]
    if w[0] == "I":
        return w[:2]
    return w[0]


def readable_profile(codes: list) -> dict:
    """The `profile` operation for DRAW profile codes, each [CODE, values...]."""
    steps, op = [], {}
    for code, *v in codes:
        if code in ("W", "WW"):
            op["outcome"] = "closed_wire" if code == "W" else "open_wire"
            break  # DRAW ignores the rest
        if code == "P":
            steps.append({"plane": {"normal": v[:3], "x": v[3:]}})
        elif code == "C":
            steps.append({"arc": {"radius": v[0], "angle": v[1]}})
        else:
            steps.append({PROFILE_STEPS[code]: v[0] if len(v) == 1 else v})
    op["steps"] = steps
    return op


class Converter:
    def __init__(self):
        self.params: dict = {}
        self.steps: list = []  # ("bind", name, op) or ("check", check)
        self.kinds: dict = {}  # name -> "shape" | "curve" | "curve2" | "plane"
        self.bop = None
        self.evolving = None

    # --- helpers ---------------------------------------------------------------------------

    def ref(self, name: str, kind: str = "shape") -> Ref:
        if name not in self.kinds:
            raise Unsupported(f"unbound {name}")
        if kind and self.kinds[name] != kind:
            raise Unsupported(f"{name} is a {self.kinds[name]}, not a {kind}")
        return Ref(name)

    def bind(self, name: str, op: dict, kind: str = "shape"):
        self.steps.append(("bind", name, op))
        self.kinds[name] = kind

    def check(self, c: dict):
        self.steps.append(("check", c))

    def nums(self, args):
        return [num(a, self.params) for a in args]

    def frame_and_values(self, args):
        """An optional leading plane, then numbers."""
        if args and args[0] in self.kinds and self.kinds[args[0]] == "plane":
            return Ref(args[0]), self.nums(args[1:])
        return None, self.nums(args)

    @staticmethod
    def with_frame(op: dict, frame):
        if frame is not None:
            op["frame"] = frame
        return op

    def curve_frame(self, v: list) -> dict:
        if len(v) == 3:
            return {"center": v}
        if len(v) == 6:
            return {"center": v[:3], "normal": v[3:]}
        if len(v) == 9:
            return {"center": v[:3], "normal": v[3:6], "x": v[6:]}
        raise Unsupported("curve frame")

    # --- commands --------------------------------------------------------------------------

    def run(self, cmd: str, a: list):
        f = getattr(self, "c_" + cmd, None)
        if f is None:
            raise Unsupported(f"command {cmd}")
        try:
            f(a)
        except (IndexError, KeyError, ValueError) as e:
            raise Unsupported(f"{cmd} arguments ({e})")

    def c_dset(self, a):
        if len(a) % 2:
            raise Unsupported("dset")
        for name, value in zip(a[::2], a[1::2]):
            v = num(value, self.params)
            if name in self.params and self.params[name] != v:
                raise Unsupported("dset again")
            self.params[name] = v

    def c_fsameparameter(self, a):
        pass

    def c_box(self, a):
        v = self.nums(a[1:])
        if len(v) == 3:
            self.bind(a[0], {"box": {"size": v}})
        elif len(v) == 6:
            self.bind(a[0], {"box": {"corner": v[:3], "size": v[3:]}})
        else:
            raise Unsupported("box form")

    def c_pcylinder(self, a):
        frame, v = self.frame_and_values(a[1:])
        op = {"radius": v[0], "height": v[1]}
        if len(v) == 3:
            op["angle"] = v[2]
        elif len(v) != 2:
            raise Unsupported("pcylinder form")
        self.bind(a[0], {"cylinder": self.with_frame(op, frame)})

    def c_pcone(self, a):
        frame, v = self.frame_and_values(a[1:])
        op = {"r1": v[0], "r2": v[1], "height": v[2]}
        if len(v) == 4:
            op["angle"] = v[3]
        elif len(v) != 3:
            raise Unsupported("pcone form")
        self.bind(a[0], {"cone": self.with_frame(op, frame)})

    def c_psphere(self, a):
        frame, v = self.frame_and_values(a[1:])
        op = {"radius": v[0]}
        if len(v) == 2:
            op["angle"] = v[1]
        elif len(v) >= 3:
            op["latitudes"] = v[1:3]
            if len(v) == 4:
                op["angle"] = v[3]
            elif len(v) != 3:
                raise Unsupported("psphere form")
        self.bind(a[0], {"sphere": self.with_frame(op, frame)})

    def c_ptorus(self, a):
        frame, v = self.frame_and_values(a[1:])
        op = {"major": v[0], "minor": v[1]}
        if len(v) == 3:
            op["angle"] = v[2]
        elif len(v) >= 4:
            op["tube"] = v[2:4]
            if len(v) == 5:
                op["angle"] = v[4]
            elif len(v) != 4:
                raise Unsupported("ptorus form")
        self.bind(a[0], {"torus": self.with_frame(op, frame)})

    def c_wedge(self, a):
        v = self.nums(a[1:])
        frame = None
        if len(v) in (13, 16):
            frame = {"origin": v[:3], "z": v[3:6], "x": v[6:9]}
            v = v[9:]
        if len(v) == 4:
            op = {"size": v[:3], "top_dx": v[3]}
        elif len(v) == 7:
            op = {"size": v[:3], "top": v[3:]}
        else:
            raise Unsupported("wedge form")
        self.bind(a[0], {"wedge": self.with_frame(op, frame)})

    def c_plane(self, a):
        v = self.nums(a[1:])
        frame = {}
        if len(v) >= 3:
            frame["origin"] = v[:3]
        if len(v) >= 6:
            frame["z"] = v[3:6]
        if len(v) == 9:
            frame["x"] = v[6:9]
        if len(v) not in (0, 3, 6, 9):
            raise Unsupported("plane form")
        self.bind(a[0], {"plane": {"frame": frame}}, "plane")

    def c_ttranslate(self, a):
        v = self.nums(a[-3:])
        for n in a[:-3]:
            self.bind(n, {"translate": {"shape": self.ref(n), "by": v}})

    def c_trotate(self, a):
        v = self.nums(a[-7:])
        for n in a[:-7]:
            self.bind(n, {"rotate": {"shape": self.ref(n), "axis": {"origin": v[:3], "direction": v[3:6]}, "angle": v[6]}})

    def c_tscale(self, a):
        v = self.nums(a[-4:])
        for n in a[:-4]:
            self.bind(n, {"scale": {"shape": self.ref(n), "center": v[:3], "factor": v[3]}})

    def c_treverse(self, a):
        for n in a:
            self.bind(n, {"reverse": {"shape": self.ref(n)}})

    def c_tcopy(self, a):
        self.bind(a[1], {"copy": {"shape": self.ref(a[0])}})

    def c_orientation(self, a):
        orient = {"F": "forward", "R": "reversed", "I": "internal", "E": "external"}.get(a[-1].upper())
        if orient is None:
            raise Unsupported("orientation")
        for n in a[:-1]:
            self.bind(n, {"orient": {"shape": self.ref(n), "as": orient}})

    def c_explode(self, a):
        pass  # the names a script uses are bound by the `nearest` calls that follow

    c_nexplode = c_explode

    def c_nearest(self, a):
        name, of, kind = a[0], a[1], a[2]
        key_of = {"E": "edge", "F": "face", "V": "vertex", "W": "wire", "SH": "shell", "SO": "solid"}
        measure_of = {"E": "length", "W": "length", "F": "area", "SH": "area", "SO": "volume"}
        sel = {"of": self.ref(of), key_of[kind]: self.nums(a[3:6])}
        if len(a) == 7 and kind in measure_of:
            sel[measure_of[kind]] = num(a[6], self.params)
        self.bind(name, {"select": sel})

    def c_polyline(self, a):
        v = self.nums(a[1:])
        self.bind(a[0], {"polyline": {"points": [v[i:i + 3] for i in range(0, len(v), 3)]}})

    def conic(self, a, op, nradii):
        v = self.nums(a[1:])
        frame, radii = v[:-nradii], v[-nradii:]
        names = {"circle": ["radius"], "ellipse": ["major", "minor"], "hyperbola": ["major", "minor"],
                 "parabola": ["focal"]}[op]
        if len(frame) in (2, 4):
            spec = {"center": frame[:2]}
            if len(frame) == 4:
                spec["x"] = frame[2:]
            kind = "curve2"
        else:
            spec = self.curve_frame(frame)
            kind = "curve"
        spec.update(dict(zip(names, radii)))
        self.bind(a[0], {op: spec}, kind)

    def c_circle(self, a):
        self.conic(a, "circle", 1)

    def c_ellipse(self, a):
        self.conic(a, "ellipse", 2)

    def c_parabola(self, a):
        self.conic(a, "parabola", 1)

    def c_beziercurve(self, a):
        n = int(a[1])
        v = self.nums(a[2:])
        if len(v) == 3 * n:
            self.bind(a[0], {"bezier": {"poles": [v[i:i + 3] for i in range(0, len(v), 3)]}}, "curve")
        elif len(v) == 4 * n:
            self.bind(a[0], {"bezier": {"poles": [v[i:i + 3] for i in range(0, len(v), 4)],
                                        "weights": v[3::4]}}, "curve")
        else:
            raise Unsupported("beziercurve form")

    def c_bsplinecurve(self, a):
        degree, nk = int(a[1]), int(a[2])
        v = self.nums(a[3:])
        knots = [[v[2 * i], v[2 * i + 1]] for i in range(nk)]
        rest = v[2 * nk:]
        if len(rest) % 4:
            raise Unsupported("bsplinecurve form")
        op = {"degree": degree, "knots": knots, "poles": [rest[i:i + 3] for i in range(0, len(rest), 4)]}
        weights = rest[3::4]
        if any(w != 1 for w in weights):
            op["weights"] = weights
        self.bind(a[0], {"bspline": op}, "curve")

    def c_mkedge(self, a):
        op = {"curve": self.ref(a[1], "curve")}
        if len(a) == 4:
            op["range"] = self.nums(a[2:])
        elif len(a) != 2:
            raise Unsupported("mkedge form")
        self.bind(a[0], {"edge": op})

    def c_vertex(self, a):
        if len(a) != 4:
            raise Unsupported("vertex form")
        self.bind(a[0], {"vertex": {"at": self.nums(a[1:])}})

    def c_edge(self, a):
        if len(a) != 3:
            raise Unsupported("edge form")
        self.bind(a[0], {"edge": {"from": self.ref(a[1]), "to": self.ref(a[2])}})

    def c_wire(self, a):
        if any(x.startswith("-") for x in a):
            raise Unsupported("wire options")
        self.bind(a[0], {"wire": {"edges": [self.ref(n) for n in a[1:]]}})

    def c_mkplane(self, a):
        if len(a) != 2:
            raise Unsupported("mkplane form")
        self.bind(a[0], {"face": {"wire": self.ref(a[1])}})

    def c_prism(self, a):
        if len(a) != 5:
            raise Unsupported("prism form")
        self.bind(a[0], {"extrude": {"profile": self.ref(a[1]), "by": self.nums(a[2:5])}})

    def c_revol(self, a):
        if len(a) != 9:
            raise Unsupported("revol form")
        v = self.nums(a[2:])
        self.bind(a[0], {"revolve": {"profile": self.ref(a[1]), "axis": {"origin": v[:3], "direction": v[3:6]},
                                     "angle": v[6]}})

    def c_pipe(self, a):
        if len(a) != 3:
            raise Unsupported("pipe form")
        self.bind(a[0], {"pipe": {"spine": self.ref(a[1]), "profile": self.ref(a[2])}})

    def c_nurbsconvert(self, a):
        if len(a) != 2:
            raise Unsupported("nurbsconvert form")
        self.bind(a[0], {"to_nurbs": {"shape": self.ref(a[1])}})

    def c_compound(self, a):
        self.bind(a[-1], {"compound": {"of": [self.ref(n) for n in a[:-1]]}})

    def c_shape(self, a):
        kind = {"SO": "solid", "SH": "shell", "C": "compound", "CS": "compsolid", "W": "wire",
                "F": "face"}.get(a[1].upper())
        if kind is None:
            raise Unsupported("shape kind")
        self.bind(a[0], {"empty": kind})

    def c_add(self, a):
        self.bind(a[1], {"add": {"to": self.ref(a[1]), "shape": self.ref(a[0])}})

    def c_halfspace(self, a):
        self.bind(a[0], {"half_space": {"face": self.ref(a[1]), "toward": self.nums(a[2:5])}})

    def c_restore(self, a):
        if len(a) != 2 or not a[0].startswith("@data:"):
            raise Unsupported("restore form")
        self.bind(a[1], {"load": {"brep": a[0][len("@data:"):]}})

    def boolean(self, op, result, obj, tool):
        self.bind(result, {op: {"object": self.ref(obj), "tool": self.ref(tool)}})

    def c_bfuse(self, a):
        self.boolean("union", *a[:3])

    def c_bcut(self, a):
        self.boolean("difference", *a[:3])

    def c_bcommon(self, a):
        self.boolean("intersection", *a[:3])

    def c_bop(self, a):
        self.ref(a[0]), self.ref(a[1])
        self.bop = (a[0], a[1])

    def bop_result(self, op, a, swap=False):
        if self.bop is None:
            raise Unsupported("no bop")
        o, t = self.bop
        self.boolean(op, a[0], *((t, o) if swap else (o, t)))

    def c_bopfuse(self, a):
        self.bop_result("union", a)

    def c_bopcut(self, a):
        self.bop_result("difference", a)

    def c_bopcommon(self, a):
        self.bop_result("intersection", a)

    def c_boptuc(self, a):
        self.bop_result("difference", a, swap=True)

    def c_bfuseblend(self, a, op="filleted_union"):
        if len(a) not in (4, 5):
            raise Unsupported("bfuseblend form")
        self.bind(a[0], {op: {"object": self.ref(a[1]), "tool": self.ref(a[2]), "radius": num(a[3], self.params)}})

    def c_bcutblend(self, a):
        self.c_bfuseblend(a, "filleted_difference")

    def c_blend(self, a):
        if len(a) < 4 or len(a) % 2:
            raise Unsupported("blend form")
        edges = [{"edge": self.ref(e), "radius": num(r, self.params)} for r, e in zip(a[2::2], a[3::2])]
        self.bind(a[0], {"fillet": {"shape": self.ref(a[1]), "edges": edges}})

    def c_mkevol(self, a):
        self.evolving = (a[0], self.ref(a[1]), [])

    def c_updatevol(self, a):
        if self.evolving is None or len(a) < 3 or (len(a) - 1) % 2:
            raise Unsupported("updatevol")
        v = self.nums(a[1:])
        self.evolving[2].append({"edge": self.ref(a[0]), "law": [v[i:i + 2] for i in range(0, len(v), 2)]})

    def c_buildevol(self, a):
        if self.evolving is None:
            raise Unsupported("buildevol")
        name, shape, edges = self.evolving
        self.bind(name, {"fillet": {"shape": shape, "edges": edges}})
        self.evolving = None

    def c_tclean(self, a):
        pass  # removes meshes

    def c_copy(self, a):
        if len(a) % 2:
            raise Unsupported("copy form")
        for src, dst in zip(a[::2], a[1::2]):
            if src not in self.kinds:
                raise Unsupported(f"unbound {src}")
            self.bind(dst, Ref(src), self.kinds[src])

    def c_renamevar(self, a):
        self.c_copy(a)

    def c_hyperbola(self, a):
        self.conic(a, "hyperbola", 2)

    def c_line(self, a):
        v = self.nums(a[1:])
        if len(v) == 4:
            self.bind(a[0], {"line": {"origin": v[:2], "direction": v[2:]}}, "curve2")
        elif len(v) == 6:
            self.bind(a[0], {"line": {"origin": v[:3], "direction": v[3:]}}, "curve")
        else:
            raise Unsupported("line form")

    def c_bsection(self, a):
        if len(a) < 3 or any(x not in ("-n2d", "-n2d1", "-n2d2", "-na", "-a") for x in a[3:]):
            raise Unsupported("bsection form")
        self.boolean("section", *a[:3])

    def c_bopsection(self, a):
        self.bop_result("section", a)

    def c_checksection(self, a):
        name = self.ref(a[0])
        if len(a) == 1:
            self.check({"exists": name})  # only prints the number of lone vertices
        elif len(a) == 3 and a[1] == "-r":
            self.check({"free_vertices": name, "eq": int(a[2])})
        else:
            raise Unsupported("checksection form")

    def c_checkfaults(self, a):
        if len(a) > 2 and a[2] == "-1":
            self.check({"exists": self.ref(a[0])})  # OCCT: unstable
        else:
            self.check({"faults": self.ref(a[0]), "source": self.ref(a[1])})

    def c_checkfreebounds(self, a):
        name, ref, kind, i = self.ref(a[0]), a[1], "closed", 2
        while i < len(a):
            if a[i] == "-type" and i + 1 < len(a) and a[i + 1] in ("closed", "opened"):
                kind = a[i + 1]
            elif a[i] == "-tol" and i + 1 < len(a) and float(a[i + 1]) <= 0:
                pass  # no sewing of the free edges
            else:
                raise Unsupported(f"checkfreebounds {a[i]}")
            i += 2
        if ref == "-1":
            self.check({"exists": name})  # OCCT: unstable
        else:
            self.check({"free_edges": name, "bounds": "closed" if kind == "closed" else "open",
                        "eq": int(ref)})

    def c_mkoffset(self, a):
        # `mkoffset result shape count step [a|i [alt]]`: offsets at step, 2 step, ... named
        # result_1, result_2, ...
        if len(a) < 4 or (len(a) > 4 and a[4] != "a"):
            raise Unsupported("mkoffset form")
        count, step = int(num(a[2], self.params)), num(a[3], self.params)
        shape = self.ref(a[1])
        for k in range(1, count + 1):
            self.bind(f"{a[0]}_{k}", {"offset_wire": {"shape": shape, "distance": step * k}})

    def c_mkoffset_fraction(self, a):
        # The offset group's end script: `mkoffset result s 1 [expr length(s) * fraction]`.
        self.bind(a[0], {"offset_wire": {"shape": self.ref(a[1]), "length_fraction": num(a[2], self.params)}})

    def c_offset_counts(self, a):
        # The offset group's end: vertices, edges and wires of the result (none of them if
        # it is to be empty, whatever the wire count says).
        v, e, w = (int(num(x, self.params)) for x in a[1:4])
        self.check({"count": self.ref(a[0]), "eq": {"vertices": v, "edges": e, "wires": w if v else 0}})

    # Offsets of shapes: parameters, then the shape (and faces left out), then the offset.
    def c_offsetparameter(self, a):
        # `offsetparameter tol c|p a|i [r|u]`: complete intersection (`c`) and arc or
        # intersection joins; `r` removes internal edges.
        if len(a) < 3 or a[1] != "c":
            raise Unsupported("offsetparameter form")
        self.offset_join = {"a": "arc", "i": "intersection"}[a[2]]
        self.offset_merge = len(a) > 3 and a[3] == "r"

    def c_offsetload(self, a):
        if len(a) < 2:
            raise Unsupported("offsetload form")
        self.offset_load = (self.ref(a[0]), num(a[1], self.params), [self.ref(f) for f in a[2:]])

    def c_offsetonface(self, a):
        raise Unsupported("offsetonface")

    def c_offsetperform(self, a):
        shape, distance, faces = self.offset_load
        join = getattr(self, "offset_join", "arc")
        if faces:
            op = {"thick_solid": {"shape": shape, "faces": faces, "thickness": distance, "join": join}}
        else:
            op = {"offset_shape": {"shape": shape, "distance": distance, "join": join}}
        if getattr(self, "offset_merge", False):
            raise Unsupported("offsetparameter r")
        self.bind(a[0], op)

    def c_thrusections(self, a):
        args = [x for x in a if x != "-safe"]
        # `-N`: the sections as they are, without OCCT's compatibility checks.
        as_given = bool(args) and args[0] == "-N"
        if as_given:
            args = args[1:]
        if len(args) < 5:
            raise Unsupported("thrusections form")
        result, solid, ruled = args[0], args[1], args[2]
        loft = {"sections": [self.ref(n) for n in args[3:]], "solid": solid == "1", "ruled": ruled == "1"}
        if as_given:
            loft["matching"] = "as_given"
        self.bind(result, {"loft": loft})

    # The general fuse builder: arguments collected, then an operation on them.
    def builder(self):
        if not hasattr(self, "objects"):
            self.objects, self.tools = [], []

    def c_bclearobjects(self, a):
        self.builder()
        self.objects = []

    def c_bcleartools(self, a):
        self.builder()
        self.tools = []

    def c_baddobjects(self, a):
        self.builder()
        self.objects += [self.ref(n) for n in a]

    def c_baddtools(self, a):
        self.builder()
        self.tools += [self.ref(n) for n in a]

    def c_bfillds(self, a):
        self.builder()

    def c_bnondestructive(self, a):
        pass

    def c_brunparallel(self, a):
        pass

    def c_bfuzzyvalue(self, a):
        if num(a[0], self.params) != 0:
            raise Unsupported("fuzzy booleans")

    def c_bbop(self, a, api=False):
        self.builder()
        op = {"0": "intersection", "1": "union", "2": "difference", "3": "difference21", "4": "section"}.get(a[1])
        if op is None or not self.objects or not self.tools:
            raise Unsupported("bbop form")
        obj, tool = (self.tools, self.objects) if op == "difference21" else (self.objects, self.tools)
        op = "difference" if op == "difference21" else op
        one = lambda l: l[0] if len(l) == 1 else l
        self.bind(a[0], {op: {"object": one(obj), "tool": one(tool)}})

    def c_bapibop(self, a):
        self.c_bbop(a, api=True)

    def c_bbuild(self, a):
        self.builder()
        if not self.objects:
            raise Unsupported("bbuild form")
        self.bind(a[0], {"partition": {"of": self.objects + self.tools}})

    def c_bapibuild(self, a):
        self.c_bbuild(a)

    def c_bsplit(self, a):
        self.builder()
        self.bind(a[0], {"split": {"objects": self.objects, "tools": self.tools}})

    def c_profile(self, a):
        codes, i = [], 1
        while i < len(a):
            cmd = profile_code(a[i])
            if cmd not in PROFILE_ARGS:
                raise Unsupported(f"profile step {cmd}")
            k = PROFILE_ARGS[cmd]
            values = a[i + 1:i + 1 + k]
            if len(values) != k:
                raise Unsupported("profile arguments")
            codes.append([cmd] + ([self.ref(values[0])] if cmd == "S" else self.nums(values)))
            i += 1 + k
        self.bind(a[0], {"profile": readable_profile(codes)})

    SEW_SWITCHES = {"s": "sewing", "a": "analysis", "c": "cutting", "n": "non_manifold", "p": "same_parameter",
                    "e": "floating_edges", "f": "faces"}

    def c_sewing(self, a):
        op, shapes, i = {}, [], 1
        while i < len(a):
            x = a[i]
            if x[:1] in "-+" and not NUMBER.match(x):
                low = x.lower()
                if low.startswith(("-mi", "+mi")):
                    op["min_tolerance"] = num(a[i + 1], self.params)
                    i += 2
                    continue
                if low.startswith(("-ma", "+ma")):
                    op["max_tolerance"] = num(a[i + 1], self.params)
                    i += 2
                    continue
                key = self.SEW_SWITCHES.get(low[1:2])
                if key is None:
                    raise Unsupported(f"sewing {x}")
                op[key] = x[0] == "+"
            elif NUMBER.match(x) or (x not in self.kinds and not shapes):
                if shapes or "tolerance" in op:
                    raise Unsupported("sewing numbers")
                op["tolerance"] = num(x, self.params)
            else:
                shapes.append(self.ref(x))
            i += 1
        if not shapes:
            raise Unsupported("sewing without shapes")
        self.bind(a[0], {"sew": {"shapes": shapes, **op}})

    def c_checkmaxtol(self, a):
        # Tolerances are the kernel's own: OCCT's value (within checkreal's 1e-4 + 1%) becomes an
        # upper bound, which a more precise kernel meets too.
        name, i, found = self.ref(a[0]), 1, False
        while i < len(a):
            opt = a[i]
            if opt == "-ref" and i + 1 < len(a):
                v = num(a[i + 1], self.params)
                if isinstance(v, Expr):
                    raise Unsupported("checkmaxtol expression")
                self.check({"tolerance": name, "le": v + 1e-4 + 0.01 * abs(v)})
            elif opt == "-min_tol" and i + 1 < len(a):
                self.check({"tolerance": name, "le": num(a[i + 1], self.params)})
            else:
                raise Unsupported(f"checkmaxtol {opt}")
            found = True
            i += 2
        if not found:
            self.check({"exists": name})

    def c_checkshape(self, a):
        if len(a) != 1:
            raise Unsupported("checkshape options")
        self.check({"valid": self.ref(a[0])})

    def c_require_var(self, a):
        self.check({"exists": Ref(a[0])})

    def c_checkprops(self, a):
        name = self.ref(a[0])
        checks, rel = [], 1e-2
        i = 1
        while i < len(a):
            opt, val = a[i], a[i + 1] if i + 1 < len(a) else None
            if val is None:
                raise Unsupported("checkprops form")
            if opt in ("-l", "-s", "-v"):
                checks.append((opt, val))
            elif opt == "-deps":
                rel = num(val, self.params)
            elif opt != "-eps":
                raise Unsupported(f"checkprops {opt}")
            i += 2
        for opt, val in checks:
            quantity = {"-l": "length", "-s": "area", "-v": "volume"}[opt]
            if val == "empty":
                self.check({quantity: name, "eq": 0})
            else:
                v = num(val, self.params)
                self.check({quantity: name, "eq": v, "rel": rel} if v != 0 else {quantity: name, "eq": 0})

    def c_checknbshapes(self, a):
        name = self.ref(a[0])
        counts = {}
        keys = {"-vertex": "vertices", "-edge": "edges", "-wire": "wires", "-face": "faces", "-shell": "shells",
                "-solid": "solids", "-compsolid": "compsolids", "-compound": "compounds", "-shape": "shapes"}
        placements = False
        i = 1
        while i < len(a):
            opt = a[i]
            if opt == "-t":
                placements = True
                i += 1
            elif opt == "-m":
                i += 2  # a message for the log
            elif opt in keys and i + 1 < len(a):
                counts[keys[opt]] = int(a[i + 1])
                i += 2
            else:
                raise Unsupported(f"checknbshapes {opt}")
        if not counts:
            raise Unsupported("checknbshapes without counts")
        c = {"count": name, "eq": counts}
        if placements:
            c["placements"] = True
        self.check(c)


def walk(v):
    """Every value nested in steps, dicts and lists."""
    if isinstance(v, dict):
        for x in v.values():
            yield from walk(x)
    elif isinstance(v, (list, tuple)):
        for x in v:
            yield from walk(x)
    else:
        yield v


def to_yaml(case_name: str, file: str, git_hash: str, body: list, notes: list, skip: str | None, timeout: float | None = None,
            expect_failure: bool = False, required_faulty: bool = False) -> str:
    """The YAML document for a translated script, or raises Unsupported."""
    conv = Converter()
    for cmd, args in parse_calls(body):
        conv.run(cmd, args)
    if conv.evolving is not None:
        raise Unsupported("mkevol without buildevol")
    # Checks after the last binding go to `checks`.
    last_bind = max((i for i, s in enumerate(conv.steps) if s[0] == "bind"), default=-1)
    build, checks = conv.steps[:last_bind + 1], [s[1] for s in conv.steps[last_bind + 1:]]
    out = [f"case: {key(case_name)}",
           "ported_from:",
           "  project: OCCT",
           f"  repo: {OCCT_REPO}",
           f"  git_hash: {git_hash}",
           f"  file: {file}"]
    if notes:
        out.append("notes:")
        out += [f"  - {json.dumps(n)}" for n in notes]
    if skip:
        out.append(f"skip: {json.dumps(skip)}")
    if timeout is not None:
        out.append(f"timeout: {timeout}")
    if expect_failure:
        out.append("expect: not_done")
    if required_faulty:
        # (The script requires a faulty shape: its `valid` checks are to fail once at least.)
        out.append("faulty: required")
    # The scripts' parameters (DRAW variables such as the blend group's SCALE1) are inlined:
    # arithmetic on numbers becomes its value, other expressions (with pi or functions) keep
    # their text with the values substituted.
    build = [inline(s, conv.params) for s in build]
    checks = [inline(c, conv.params) for c in checks]
    out.append("build:")
    for s in build:
        if s[0] == "bind":
            out.append(f"  - {key(s[1])}: {flow(s[2])}")
        else:
            out.append(f"  - check: {flow(s[1])}")
    if checks:
        out.append("checks:")
        out += [f"  - {flow(c)}" for c in checks]
    return "\n".join(out) + "\n"
