#!/usr/bin/env python3
"""Convert simple OCCT DRAW (Tcl) test scripts into Rust tests.

OCCT's tests/ tree holds DRAW scripts: tests/<group>/<grid>/<case>, plus begin/end files
at the group and grid level. This tool converts the *straight-line* ones (plain command
sequences, `set` with constant values, `[locate_data_file ...]`, constant `[expr ...]`, and
nested command substitution) into Rust tests. Scripts using real Tcl control flow (if, foreach,
proc, regexp, catch, ...) are skipped and listed in the report.

Each DRAW command becomes a method on `occt_tests::draw::Draw` taking its arguments as strings,
mirroring how DRAW itself dispatches commands. Missing command methods are appended to
crates/occt-tests/src/draw/commands.rs as shims that return `DrawError::NotImplemented`;
existing methods (i.e. ones you have implemented) are never modified.

Usage:
    tools/draw2rs.py [--occt ~/OCCT] [--data-dir DIR ...] [--groups boolean,offset,...]

Re-running is safe: generated test files are rewritten, shims are only appended.

Intentional differences from OCCT live in crates/occt-tests/draw_deviations.json. An entry
either replaces generated statements of a script (the reason is emitted as a comment on the
test; a stale entry, whose statement is no longer generated, aborts the run) or skips the
script altogether (the reason shows in the report's skip reasons).
"""

from __future__ import annotations

import argparse
import json

import drawyaml
import subprocess
import ast
import math
import os
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
TESTS_CRATE = REPO / "crates" / "occt-tests"
COMMANDS_RS = TESTS_CRATE / "src" / "draw" / "commands.rs"

# --- scope ---------------------------------------------------------------------------------

# Groups covering computational geometry and topology. Excluded: visualization (v3d, vselect,
# opengl*, vtk, hlr), data exchange (de*, xcaf, gdt, xml), OCAF (caf, persist, metadata),
# demos and performance tests.
DEFAULT_GROUPS = [
    "blend", "boolean", "bugs", "chamfer", "cr", "draft", "evolved", "feat", "fillet2d",
    "geometry", "heal", "helix", "lowalgos", "mesh", "mkface", "nproject", "offset", "pipe",
    "prism", "prj", "sewing", "thrusection",
]
EXCLUDED_GRIDS = {
    "bugs": {"caf", "demo", "iges", "step", "stlvrml", "xde"},
}

# Commands that only affect display, timing, logging or the DRAW session itself: dropped.
DROPPED_COMMANDS = {
    "puts", "pload", "cpulimit", "checkview", "smallview", "fit", "2dfit", "donly", "display",
    "erase", "clear", "axo", "top", "front", "left", "right", "bottom", "back", "xwd", "don",
    "isos", "protect", "dchrono", "chrono", "dlog", "decho", "dtracelevel", "view", "wait",
    "meminfo", "checkcolor", "trinfo_display", "dversion", "dsetsignal", "set_display",
    "autodisplay", "repaint", "dflush", "dbreak", "datadir", "checkplatform", "vinit",
    # (Panning and magnifying the 3D and 2D views.)
    "pu", "pd", "pl", "pr", "mu", "md", "mu4", "md4", "2dpu", "2dpd", "2dpl", "2dpr", "2dmu",
    "2dmd", "2dzoom", "zoom", "wzoom",
    # (Drawing normals of faces.)
    "normals",
    "vdisplay", "vfit", "vsetdispmode", "vdefaults", "vclear", "vviewparams", "vtop", "vaxo",
    "vbottom", "vzfit", "vright", "vfront", "vaspects", "vleft", "vclose", "vback", "vsetcolor",
    "vrotate", "vrenderparams", "verase", "vsettransparency", "vzoom", "vpan", "vsetmaterial",
    "veraseall", "vglinfo", "vcamera", "vori", "vbounding", "vcolorscale", "vtexture",
    "vtexrepeat", "vselprops", "vselmode", "vtriangle", "vsegment", "vplane", "vlocation",
    "vcircle", "vpoint", "vdrawsph", "vdrawparray", "dmode", "triedron", "nbiso", "discretisation",
    "chkfi", "tclean_display", "checktrinfo_display", "dtyp", "whatis", "dump", "pick",
    "checkdisplay",
}
# Commands whose result depends on an interactive viewer: the whole script is skipped.
VIEWER_COMMANDS = {
    "vselect", "vmoveto", "vreadpixel", "vdump", "vhlr", "vhlrtype", "vcomputehlr", "vstate",
    "vdetected", "vchangeselected", "vmanipulator", "vprs", "vdisplayall",
}
# Data-exchange commands (STEP/IGES/STL/VRML/XDE): out of scope, the script is skipped.
DE_COMMANDS = {
    "stepread", "stepwrite", "testreadstep", "testwritestep", "ReadStep", "WriteStep",
    "igesbrep", "igesread", "igeswrite", "brepiges", "testreadiges", "ReadIges", "WriteIges",
    "readstl", "writestl", "loadvrml", "writevrml", "ReadGltf", "WriteGltf", "ReadObj",
    "WriteObj", "ReadPly", "WritePly", "XNewDoc", "XOpen", "XSave", "XShow", "XGetOneShape",
    "newmodel", "tpcompound", "tpstat", "data", "param", "xstep", "XReadStep", "XWriteStep",
    "readfile", "writefile",
}
# Tcl built-ins (other than `set`, `puts` and constant `expr`) mean the script is not simple.
TCL_BUILTINS = {
    "if", "else", "elseif", "for", "foreach", "while", "proc", "catch", "regexp", "regsub",
    "lindex", "llength", "lappend", "lrange", "lsearch", "lsort", "lreplace", "linsert", "lset",
    "list", "string", "format", "scan", "incr", "append", "array", "info", "global", "upvar",
    "uplevel", "return", "break", "continue", "switch", "error", "eval", "source", "file",
    "open", "close", "gets", "read", "exec", "split", "join", "concat", "unset", "after",
    "exit", "binary", "clock", "namespace", "subst", "dict", "encoding", "glob", "pwd", "cd",
    "variable", "update", "trace", "rename", "interp", "package", "fconfigure", "flush",
    "seek", "tell", "eof", "socket", "vwait", "time", "expr", "lassign", "lreverse", "lmap",
    "apply", "tailcall", "try", "throw", "chan", "zlib", "fcopy", "load",
}

# Names reserved by the hand-written part of `Draw`; command methods get a `_cmd` suffix.
RESERVED_METHODS = {"new", "get", "set", "data_file", "not_implemented", "default", "require_var",
                    "shape_var", "surface_var", "set_shape_var", "model"}

RUST_KEYWORDS = {
    "as", "break", "const", "continue", "crate", "else", "enum", "extern", "false", "fn", "for",
    "if", "impl", "in", "let", "loop", "match", "mod", "move", "mut", "pub", "ref", "return",
    "self", "Self", "static", "struct", "super", "trait", "true", "type", "unsafe", "use",
    "where", "while", "async", "await", "dyn", "abstract", "become", "box", "do", "final",
    "macro", "override", "priv", "typeof", "unsized", "virtual", "yield", "try", "gen", "union",
}


class Unsupported(Exception):
    """The script uses something outside the simple subset; message is a short reason."""


# --- Tcl lexing ----------------------------------------------------------------------------

@dataclass
class Lit:
    text: str


@dataclass
class Var:
    name: str


@dataclass
class Sub:
    """Command substitution `[...]`."""
    commands: list  # list[list[Word]]


@dataclass
class Word:
    parts: list  # list[Lit | Var | Sub]
    braced: bool = False

    def literal(self) -> str | None:
        if all(isinstance(p, Lit) for p in self.parts):
            return "".join(p.text for p in self.parts)
        return None


@dataclass
class Command:
    words: list  # list[Word]
    line: int


_BACKSLASH = {"n": "\n", "t": "\t", "r": "\r", "a": "\a", "b": "\b", "f": "\f", "v": "\v"}


class TclLexer:
    """Tcl's dodekalogue, minus `{*}` (rejected) and with `$a(idx)` arrays rejected."""

    def __init__(self, src: str):
        self.s = src
        self.i = 0
        self.line = 1

    def peek(self, k: int = 0) -> str:
        j = self.i + k
        return self.s[j] if j < len(self.s) else ""

    def adv(self) -> str:
        ch = self.s[self.i]
        self.i += 1
        if ch == "\n":
            self.line += 1
        return ch

    def script(self, close: str = "") -> list[Command]:
        cmds: list[Command] = []
        while True:
            # skip whitespace and command separators
            while self.i < len(self.s):
                ch = self.peek()
                if ch in " \t\r\n;":
                    self.adv()
                elif ch == "\\" and self.peek(1) == "\n":
                    self.adv()
                    self.adv()
                else:
                    break
            if self.i >= len(self.s):
                if close:
                    raise Unsupported("unterminated [")
                return cmds
            if close and self.peek() == close:
                self.adv()
                return cmds
            if self.peek() == "#":
                self.skip_comment()
                continue
            cmds.append(self.command(close))

    def skip_comment(self) -> None:
        while self.i < len(self.s):
            ch = self.adv()
            if ch == "\\" and self.i < len(self.s):
                self.adv()
            elif ch == "\n":
                return

    def command(self, close: str) -> Command:
        line = self.line
        words: list[Word] = []
        while self.i < len(self.s):
            ch = self.peek()
            if ch in " \t\r":
                self.adv()
                continue
            if ch == "\\" and self.peek(1) == "\n":
                self.adv()
                self.adv()
                continue
            if ch in "\n;" or (close and ch == close):
                break
            words.append(self.word(close))
        return Command(words, line)

    def word(self, close: str) -> Word:
        ch = self.peek()
        if ch == "{":
            if self.s.startswith("{*}", self.i):
                raise Unsupported("{*} expansion")
            return Word([Lit(self.braced())], braced=True)
        if ch == '"':
            self.adv()
            parts = self.parts(lambda c: c == '"', close)
            if self.peek() != '"':
                raise Unsupported("unterminated quote")
            self.adv()
            nxt = self.peek()
            if nxt and nxt not in " \t\r\n;" and not (close and nxt == close):
                raise Unsupported("extra characters after close-quote")
            return Word(parts)
        return Word(self.parts(lambda c: c in " \t\r\n;" or (bool(close) and c == close), close))

    def braced(self) -> str:
        assert self.adv() == "{"
        depth = 1
        out = []
        while self.i < len(self.s):
            ch = self.adv()
            if ch == "\\" and self.i < len(self.s):
                nxt = self.adv()
                if nxt == "\n":
                    out.append(" ")
                else:
                    out.append(ch + nxt)
                continue
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    return "".join(out)
            out.append(ch)
        raise Unsupported("unterminated {")

    def parts(self, stop, close: str) -> list:
        parts: list = []
        buf: list[str] = []

        def flush():
            if buf:
                parts.append(Lit("".join(buf)))
                buf.clear()

        while self.i < len(self.s):
            ch = self.peek()
            if stop(ch):
                break
            # (A backslash-newline ends a bare word: it reads as a space between words.)
            if ch == "\\" and self.peek(1) == "\n" and stop(" "):
                break
            if ch == "\\":
                self.adv()
                nxt = self.adv() if self.i < len(self.s) else ""
                if nxt == "\n":
                    buf.append(" ")
                else:
                    buf.append(_BACKSLASH.get(nxt, nxt))
                continue
            if ch == "[":
                self.adv()
                flush()
                parts.append(Sub(self.script("]")))
                continue
            if ch == "$":
                name = self.varname()
                if name is None:
                    buf.append(self.adv())
                else:
                    flush()
                    parts.append(Var(name))
                continue
            buf.append(self.adv())
        flush()
        return parts

    def varname(self) -> str | None:
        assert self.peek() == "$"
        if self.peek(1) == "{":
            end = self.s.find("}", self.i)
            if end < 0:
                raise Unsupported("unterminated ${")
            name = self.s[self.i + 2:end]
            self.i = end + 1
            return name
        m = re.compile(r"(?:::)?[A-Za-z0-9_]+(?:::[A-Za-z0-9_]+)*").match(self.s, self.i + 1)
        if not m:
            return None
        self.i = m.end()
        if self.peek() == "(":
            raise Unsupported("array variable")
        return m.group(0)


def parse_tcl(src: str) -> list[Command]:
    return TclLexer(src).script()


# --- static evaluation of `expr` -----------------------------------------------------------

_EXPR_FUNCS = {
    "sqrt": math.sqrt, "sin": math.sin, "cos": math.cos, "tan": math.tan, "asin": math.asin,
    "acos": math.acos, "atan": math.atan, "atan2": math.atan2, "exp": math.exp, "log": math.log,
    "log10": math.log10, "pow": math.pow, "abs": abs, "fabs": abs, "double": float,
    "int": lambda x: int(x), "round": lambda x: int(math.floor(x + 0.5)), "floor": math.floor,
    "ceil": math.ceil, "hypot": math.hypot, "sinh": math.sinh, "cosh": math.cosh,
    "tanh": math.tanh, "min": min, "max": max, "fmod": math.fmod,
}


def _number(text: str):
    t = text.strip()
    try:
        if re.fullmatch(r"[+-]?\d+", t):
            return int(t)
        return float(t)
    except ValueError:
        return None


def eval_expr(text: str):
    """Evaluate a constant Tcl expression; raises Unsupported if not constant arithmetic."""
    t = text.replace("&&", " and ").replace("||", " or ")
    t = re.sub(r"(?<![=!<>])!(?!=)", " not ", t)
    try:
        tree = ast.parse(t.strip(), mode="eval")
    except SyntaxError:
        raise Unsupported("non-constant expr")

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return n.value
        if isinstance(n, ast.UnaryOp):
            v = ev(n.operand)
            if isinstance(n.op, ast.USub):
                return -v
            if isinstance(n.op, ast.UAdd):
                return v
            if isinstance(n.op, ast.Not):
                return int(not v)
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add):
                return a + b
            if isinstance(n.op, ast.Sub):
                return a - b
            if isinstance(n.op, ast.Mult):
                return a * b
            if isinstance(n.op, ast.Div):
                if isinstance(a, int) and isinstance(b, int):
                    return a // b
                return a / b
            if isinstance(n.op, ast.Mod):
                return a % b
            if isinstance(n.op, ast.Pow):
                return a ** b
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in _EXPR_FUNCS:
            return _EXPR_FUNCS[n.func.id](*[ev(a) for a in n.args])
        if isinstance(n, ast.Name) and n.id == "pi":
            return math.pi
        raise Unsupported("non-constant expr")

    try:
        return ev(tree)
    except (ZeroDivisionError, ValueError, TypeError, OverflowError):
        raise Unsupported("non-constant expr")


def tcl_number_str(v) -> str:
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, int):
        return str(v)
    if v == int(v) and abs(v) < 1e15:
        return f"{v:.1f}"
    return repr(v)


# --- translation ---------------------------------------------------------------------------

@dataclass
class RVal:
    """A Rust `&str`-valued expression."""
    code: str
    const: str | None = None  # known constant value


@dataclass
class CaseResult:
    rust_body: list = field(default_factory=list)
    commands: Counter = field(default_factory=Counter)
    examples: dict = field(default_factory=dict)
    data_files: list = field(default_factory=list)
    todos: list = field(default_factory=list)


def rust_str(s: str) -> str:
    out = s.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\t", "\\t")
    out = out.replace("\r", "\\r")
    return f'"{out}"'


def cmd_method(name: str) -> str:
    m = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if m[0].isdigit():
        m = "_" + m
    if m in RESERVED_METHODS:
        m += "_cmd"
    if m in RUST_KEYWORDS:
        m = "r#" + m
    return m


def rust_ident(name: str) -> str:
    base = re.sub(r"[^A-Za-z0-9_]", "_", name).lower() or "v"
    if base[0].isdigit():
        base = "v_" + base
    if base in RUST_KEYWORDS or base in ("d", "self", "super", "crate") or re.fullmatch(r"t\d+", base):
        base = base + "_"
    return base


class Translator:
    def __init__(self, consts: dict, procs: set, data_index: dict):
        self.consts = dict(consts)  # Tcl var -> constant string value
        self.rvars: dict[str, str] = {}  # Tcl var -> Rust local holding a String
        self.procs = procs
        self.data_index = data_index
        self.res = CaseResult()
        self.tmp = 0
        self.locals: set[str] = set()

    # value of a word, possibly hoisting nested commands into temporaries
    def word(self, w: Word, pre: list[str]) -> RVal:
        lit = w.literal()
        if lit is not None:
            return RVal(rust_str(lit), lit)
        vals = [self.part(p, pre) for p in w.parts]
        if all(v.const is not None for v in vals):
            s = "".join(v.const for v in vals)
            return RVal(rust_str(s), s)
        if len(vals) == 1:
            return vals[0]
        fmt, args = "", []
        for v in vals:
            if v.const is not None:
                fmt += v.const.replace("{", "{{").replace("}", "}}")
            else:
                fmt += "{}"
                args.append(v.code)
        return RVal(f"format!({rust_str(fmt)}, {', '.join(args)}).as_str()")

    def part(self, p, pre: list[str]) -> RVal:
        if isinstance(p, Lit):
            return RVal(rust_str(p.text), p.text)
        if isinstance(p, Var):
            if p.name in self.consts:
                v = self.consts[p.name]
                return RVal(rust_str(v), v)
            if p.name in self.rvars:
                return RVal(f"{self.rvars[p.name]}.as_str()")
            raise Unsupported("unknown variable")
        if isinstance(p, Sub):
            if len(p.commands) != 1:
                raise Unsupported("multi-command substitution")
            return self.subst(p.commands[0], pre)
        raise AssertionError(p)

    def subst(self, c: Command, pre: list[str]) -> RVal:
        name = c.words[0].literal() if c.words else None
        if name == "locate_data_file" and len(c.words) == 2:
            fname = c.words[1].literal()
            if fname is None:
                raise Unsupported("dynamic data file name")
            self.res.data_files.append(fname)
            return RVal(f"data_file({rust_str(fname)}).as_str()")
        if name == "expr":
            text = " ".join(self.word(w, pre).const or "\0" for w in c.words[1:])
            if "\0" in text:
                raise Unsupported("non-constant expr")
            v = tcl_number_str(eval_expr(text))
            return RVal(rust_str(v), v)
        call = self.call(c, pre)
        self.tmp += 1
        t = f"t{self.tmp}"
        pre.append(f"let {t} = {call}?;")
        return RVal(f"{t}.as_str()")

    def call(self, c: Command, pre: list[str]) -> str:
        name = c.words[0].literal()
        if name is None:
            raise Unsupported("dynamic command name")
        self.check_command(name)
        args = [self.word(w, pre) for w in c.words[1:]]
        self.res.commands[name] += 1
        if name not in self.res.examples:
            self.res.examples[name] = " ".join(
                (w.literal() if w.literal() is not None else "…") for w in c.words)[:100]
        return f"d.{cmd_method(name)}(&[{', '.join(a.code for a in args)}])"

    def check_command(self, name: str) -> None:
        if name in TCL_BUILTINS:
            raise Unsupported(f"tcl:{name}")
        if name in self.procs:
            raise Unsupported(f"proc:{name}")
        if name in VIEWER_COMMANDS:
            raise Unsupported(f"viewer:{name}")
        if name in DE_COMMANDS:
            raise Unsupported(f"data-exchange:{name}")
        if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", name):
            raise Unsupported("odd command name")

    def statement(self, c: Command) -> None:
        if not c.words:
            return
        name = c.words[0].literal()
        if name is None:
            raise Unsupported("dynamic command name")
        pre: list[str] = []
        if name == "puts":
            self.puts(c, pre)
        elif name == "set":
            self.set(c, pre)
            return
        elif name in DROPPED_COMMANDS:
            # still evaluate command substitutions in `puts [checkshape r]`-style lines
            return
        elif name == "repeat":
            # `repeat N body`: the body's statements in a loop (dropped if none remain).
            count = c.words[1].literal() if len(c.words) == 3 else None
            body = c.words[2].literal() if len(c.words) == 3 else None
            if count is None or body is None or not count.isdigit():
                raise Unsupported("dynamic repeat")
            outer = self.res.rust_body
            self.res.rust_body = []
            for cc in parse_tcl(body):
                self.statement(cc)
            inner, self.res.rust_body = self.res.rust_body, outer
            if inner:
                self.res.rust_body.append(f"for _ in 0..{count} {{")
                self.res.rust_body.extend("    " + line for line in inner)
                self.res.rust_body.append("}")
            return
        elif name in PROC_EXPANSIONS and name in self.procs:
            # A known proc of the test groups: its body, with the arguments substituted.
            args = [self.word(w, pre).const for w in c.words[1:]]
            if any(a is None for a in args):
                raise Unsupported(f"proc:{name} with computed arguments")
            for cc in parse_tcl(PROC_EXPANSIONS[name](*args)):
                self.statement(cc)
            return
        else:
            pre_call = self.call(c, pre)
            pre.append(f"{pre_call}?;")
        self.res.rust_body.extend(pre)

    def puts(self, c: Command, pre: list[str]) -> None:
        if len(c.words) == 2:
            lit = c.words[1].literal()
            if lit is not None:
                m = re.match(r"TODO\s+(\S+)\s+([^:]+):\s*(.*)", lit)
                if m:
                    self.res.todos.append((m.group(1), m.group(2).strip(), m.group(3)))
                return
            parts = c.words[1].parts
            if len(parts) == 1 and isinstance(parts[0], Sub):
                rv = self.subst(parts[0].commands[0], pre)
                if rv.const is None and pre and pre[-1].startswith("let t"):
                    pre[-1] = pre[-1].replace(f"let t{self.tmp} = ", "")
                    self.tmp -= 1
                return
        for w in c.words[1:]:
            for p in w.parts:
                if isinstance(p, Sub):
                    self.subst(p.commands[0], pre)

    def set(self, c: Command, pre: list[str]) -> None:
        if len(c.words) != 3:
            raise Unsupported("tcl:set (read)")
        var = c.words[1].literal()
        if var is None or "(" in var:
            raise Unsupported("dynamic set")
        val = self.word(c.words[2], pre)
        self.rvars.pop(var, None)
        self.consts.pop(var, None)
        if val.const is not None:
            self.consts[var] = val.const
            self.res.rust_body.extend(pre)
            return
        local = rust_ident(var)
        self.locals.add(local)
        tmp = f"t{self.tmp}"
        if val.code == f"{tmp}.as_str()" and pre and pre[-1].startswith(f"let {tmp} = "):
            # `set x [cmd ...]`: bind the command result directly
            pre[-1] = pre[-1].replace(f"let {tmp} = ", f"let {local} = ", 1)
            self.tmp -= 1
        elif re.fullmatch(r"(data_file|format!)\(.*\)\.as_str\(\)", val.code):
            pre.append(f"let {local} = {val.code[:-len('.as_str()')]};")
        else:
            pre.append(f"let {local} = {val.code}.to_string();")
        self.rvars[var] = local
        self.res.rust_body.extend(pre)


# --- OCCT test tree ------------------------------------------------------------------------

@dataclass
class Case:
    group: str
    grid: str
    name: str
    path: Path


def discover(occt: Path, groups: list[str]) -> list[Case]:
    cases: list[Case] = []
    tests = occt / "tests"
    for g in groups:
        gdir = tests / g
        if not gdir.is_dir():
            print(f"warning: no group {g}", file=sys.stderr)
            continue
        grids_list = gdir / "grids.list"
        grids = []
        if grids_list.exists():
            for line in grids_list.read_text().splitlines():
                parts = line.split()
                if len(parts) >= 2 and not line.startswith("#"):
                    grids.append(parts[1])
        else:
            grids = sorted(p.name for p in gdir.iterdir() if p.is_dir())
        for grid in grids:
            if grid in EXCLUDED_GRIDS.get(g, set()):
                continue
            gr = gdir / grid
            if not gr.is_dir():
                continue
            for f in sorted(gr.iterdir()):
                if f.is_file() and f.name not in ("begin", "end", "cases.list") and not f.name.startswith("."):
                    cases.append(Case(g, grid, f.name, f))
    return cases


def preamble_info(path: Path) -> tuple[dict, set, bool, list]:
    """Constants and proc names from a begin file; flag whether it does anything else that
    matters (assignments from commands etc.) - in that case we can't model it. Also its
    `dset` commands (DRAW variables), which each script must run first."""
    consts: dict = {}
    procs: set = set()
    complex_ = False
    dsets: list = []
    if not path.exists():
        return consts, procs, complex_, dsets
    try:
        cmds = parse_tcl(path.read_text(errors="replace"))
    except Unsupported:
        return consts, procs, True, dsets
    for c in cmds:
        name = c.words[0].literal() if c.words else None
        if name == "proc" and len(c.words) >= 2:
            procs.add(c.words[1].literal() or "")
        elif name == "set" and len(c.words) == 3:
            var, val = c.words[1].literal(), c.words[2].literal()
            if var and val is not None:
                consts[var] = val
        elif name == "dset":
            dsets.append(c)
        elif name in ("if", "pload", "cpulimit", "protect", "puts", "global") or name in DROPPED_COMMANDS:
            continue
        else:
            complex_ = True
    return consts, procs, complex_, dsets


def build_data_index(dirs: list[Path]) -> dict:
    idx: dict = {}
    for d in dirs:
        if not d.is_dir():
            continue
        for root, _dirs, files in os.walk(d):
            for f in files:
                idx.setdefault(f, os.path.join(root, f))
    return idx


def _offsetshape(distance: str, faces: str, calcul: str, join: str) -> str:
    """tests/offset/begin's OFFSETSHAPE: offsets the shape `s` (leaving out `faces`, of `s`)."""
    faces = faces.split()
    return ((f"explode s f\n" if faces else "")
            + f"offsetparameter 1e-7 {calcul} {join}\noffsetload s {distance} {' '.join(faces)}\noffsetperform result\n")


# Procs of the test groups expanded where called, by their bodies.
PROC_EXPANSIONS = {"OFFSETSHAPE": _offsetshape}


# Group epilogues modelled on tests/<group>/end (which run when the script defines `result`):
# most check that `result` exists and is valid. The ends of UNTRANSLATED_END do the test's
# real work (offset runs the offset algorithm and checks it, mesh meshes, ...); a case of
# theirs with no checks of its own is not ported. The bugs end only checks meshes.
_CHECK_RESULT = ['d.checkshape(&["result"])?;']
GROUP_EPILOGUE = {
    g: _CHECK_RESULT
    for g in ("boolean", "blend", "chamfer", "draft", "feat", "mkface", "nproject", "prism",
              "prj", "thrusection")
}
GROUP_EPILOGUE["sewing"] = ['d.require_var("result")?;']
GROUP_EPILOGUE["heal"] = ['d.require_var("result")?;']
UNTRANSLATED_END = {"offset", "mesh", "helix", "fillet2d"}
# Grid ends (tests/<group>/<grid>/end) run after each case of the grid; those not modelled here
# are translated like the case, and a case whose grid end does not translate is not ported.
GRID_EPILOGUE = {
    # "If result exists, check it; else the pipe cannot be built."
    f"pipe/{g}": ['d.checkshape(&["result"])?;'] for g in ("standard", "bugs", "specific")
}
# The offset group's end script, for its grids offsetting wires (`command mkoffset`): an offset
# by a fraction of the wire's length, its length and sub-shape counts checked (`offset_counts`:
# the end counts the result itself as a wire when it has none, so an empty result is expected
# to have no wires).
_WIRE_OFFSET_END = """
mkoffset_fraction result s $off_param
checkprops result -l $length
offset_counts result $nbsh_v $nbsh_e $nbsh_w
checkshape result
"""
GRID_EPILOGUE_TCL = {
    f"offset/wire_{kind}_{size}": _WIRE_OFFSET_END
    for kind in ("closed_inside", "closed_outside", "unclosed_outside") for size in ("0_005", "0_025", "0_075")
}
# The volume maker's end: `bopcheck result` must find the result OK.
GRID_EPILOGUE["boolean/volumemaker"] = ['d.require_no_interference("result")?;']
# Ends that only draw the result for OCCT's image comparison.
GRID_EPILOGUE.update({f"boolean/{g}_2d": [] for g in ("bcommon", "bcut", "bfuse", "bopcommon", "bopcut",
                                                      "bopfuse", "boptuc")})

# --- emission ------------------------------------------------------------------------------

def rust_test_name(case: str, used: set) -> str:
    n = re.sub(r"[^A-Za-z0-9_]", "_", case).lower()
    if not n or n[0].isdigit():
        n = "case_" + n
    if n in RUST_KEYWORDS:
        n = "r#" + n
    base, k = n, 2
    while n in used:
        n = f"{base}_{k}"
        k += 1
    used.add(n)
    return n


def todo_applies(platforms: str) -> bool:
    for p in re.split(r"[\s,]+", platforms):
        if p in ("ALL", "Linux", "Debian", "Ubuntu"):
            return True
    return False


def translate_case(case: Case, consts: dict, procs: set, data_index: dict, prelude: list = ()) -> CaseResult:
    src = case.path.read_text(errors="replace")
    cmds = list(prelude) + parse_tcl(src)
    tr = Translator(consts, procs, data_index)
    for c in cmds:
        tr.statement(c)
    grid_end = case.path.parent / "end"
    if f"{case.group}/{case.grid}" in GRID_EPILOGUE_TCL:
        for c in parse_tcl(GRID_EPILOGUE_TCL[f"{case.group}/{case.grid}"]):
            tr.statement(c)
    elif f"{case.group}/{case.grid}" in GRID_EPILOGUE:
        tr.res.rust_body.extend(GRID_EPILOGUE[f"{case.group}/{case.grid}"])
    elif grid_end.is_file():
        try:
            for c in parse_tcl(grid_end.read_text(errors="replace")):
                tr.statement(c)
        except Unsupported as e:
            raise Unsupported(f"grid end file: {e}")
    if case.group in GROUP_EPILOGUE and re.search(r"\bresult\b", src):
        tr.res.rust_body.extend(GROUP_EPILOGUE[case.group])
    return tr.res


DEVIATIONS_JSON = TESTS_CRATE / "draw_deviations.json"
EXPECTED_FAILURES = TESTS_CRATE / "expected_failures.txt"


def load_expected_failures() -> set:
    """Scripts whose tests are expected to fail (known gaps in kentron): their YAML cases are
    `expect: not_done` and their Rust tests `#[should_panic]`, passing while they fail."""
    if not EXPECTED_FAILURES.exists():
        return set()
    lines = (l.split("#", 1)[0].strip() for l in EXPECTED_FAILURES.read_text().splitlines())
    return {l for l in lines if l}


def load_deviations() -> dict:
    """Intentional differences from OCCT: script path -> [{from, to, reason}] (statements
    replaced), [{note}] (known, left failing), [{timeout, reason}] (a shorter processor-time
    limit for a case known to run long and fail), or {"skip": reason} (the script is not
    converted)."""
    if not DEVIATIONS_JSON.exists():
        return {}
    data = json.loads(DEVIATIONS_JSON.read_text())
    return {k: v for k, v in data.items() if not k.startswith("_")}


def apply_deviations(rel: str, body: list, deviations: dict) -> list:
    """Replaces the deviating statements of a script; returns the reasons. Fails on a
    deviation whose statement is not generated any more (a stale entry)."""
    reasons = []
    entry = deviations.get(rel, [])
    for dev in entry if isinstance(entry, list) else []:
        if "note" in dev:
            reasons.append(dev["note"])  # a known difference, the test left as it is
            continue
        if "timeout" in dev:
            reasons.append(dev["reason"])  # a shorter limit (see `case_timeout`)
            continue
        if dev["from"] not in body:
            raise SystemExit(f"stale deviation for {rel}: statement not found: {dev['from']}")
        body[body.index(dev["from"])] = dev["to"]
        reasons.append(dev["reason"])
    return reasons


SELECTORS_JSON = TESTS_CRATE / "draw_selectors.json"


def load_selectors() -> dict:
    """Sub-shapes that scripts name by OCCT's enumeration order (`explode s E` then `s_3`),
    pinned down by position: script path -> [{"explode": k, "names": {name: [type, x, y,
    z, measure]}}] (the centroid, and the length or area that tells apart sub-shapes sharing
    one), k counting the script's explode/nexplode calls from 0. Written by
    tools/occt_selectors.py from runs of OCCT's DRAW; kentron enumerates differently."""
    if not SELECTORS_JSON.exists():
        return {}
    data = json.loads(SELECTORS_JSON.read_text())
    return {k: v for k, v in data.items() if not k.startswith("_")}


EXPLODE_RE = re.compile(r'^d\.n?explode\(&\["([^"]*)"')


def apply_selectors(rel: str, body: list, selectors: dict) -> list:
    """After each explode with recorded selectors, rebinds the names the script uses to the
    sub-shapes at OCCT's positions (kentron's `nearest` command)."""
    entries = {e["explode"]: e["names"] for e in selectors.get(rel, [])}
    if not entries:
        return body
    out, k = [], 0
    for stmt in body:
        out.append(stmt)
        m = EXPLODE_RE.match(stmt)
        if m:
            for name, (kind, x, y, z, measure) in sorted(entries.get(k, {}).items()):
                out.append(f'd.nearest(&["{name}", "{m.group(1)}", "{kind}", "{x!r}", "{y!r}", "{z!r}", "{measure!r}"])?;')
            k += 1
    return out


def case_timeout(rel: str, deviations: dict):
    """The processor-time limit a deviation sets for a script's case, if any."""
    entry = deviations.get(rel, [])
    return next((dev["timeout"] for dev in entry if "timeout" in dev), None) if isinstance(entry, list) else None


def try_yaml(case: Case, res: CaseResult, data_index: dict, occt: Path, deviations: dict, selectors: dict,
             git_hash: str, expected: set):
    """The case as a YAML document (docs/test-cases.md), or None if it needs Rust."""
    rel = case.path.relative_to(occt).as_posix()
    body = list(res.rust_body)
    notes = [f"kentron deviation: {r}" for r in apply_deviations(rel, body, deviations)]
    if rel in selectors:
        body = apply_selectors(rel, body, selectors)
    notes += [f"OCCT TODO {tid} ({plat}): {msg}" for tid, plat, msg in res.todos]
    known_fail = [t for t in res.todos if todo_applies(t[1]) and not t[1].startswith("?")]
    # Missing data is decided when the cases run (OCCT's dataset is an external download).
    skip = None
    if known_fail:
        skip = f"known failure in OCCT: {known_fail[0][0]}"
    try:
        return drawyaml.to_yaml(case.name, rel, git_hash, body, notes, skip, case_timeout(rel, deviations),
                                rel in expected)
    except drawyaml.Unsupported as e:
        if os.environ.get("YAML_REASONS"):
            with open(os.environ["YAML_REASONS"], "a") as f:
                f.write(f"{str(e)[:60]}\t{rel}\n")
        return None


def emit_case(case: Case, res: CaseResult, fn_name: str, data_index: dict, occt: Path,
              deviations: dict, selectors: dict, expected: set) -> str:
    rel = case.path.relative_to(occt)
    lines = [f"/// `{rel.as_posix()}`"]
    for reason in apply_deviations(rel.as_posix(), res.rust_body, deviations):
        lines.append(f"/// kentron deviation (draw_deviations.json): {reason}")
    if rel.as_posix() in selectors:
        res.rust_body = apply_selectors(rel.as_posix(), res.rust_body, selectors)
        lines.append("/// kentron: sub-shapes named by OCCT's enumeration are selected by position (draw_selectors.json)")
    missing = sorted({f for f in res.data_files if f not in data_index})
    known_fail = [t for t in res.todos if todo_applies(t[1]) and not t[1].startswith("?")]
    for tid, plat, msg in res.todos:
        lines.append(f"/// OCCT TODO {tid} ({plat}): {msg}")
    if missing:
        lines.append(f'#[ignore = "missing OCCT test data: {", ".join(missing)}"]')
    elif known_fail:
        lines.append(f'#[ignore = "known failure in OCCT: {known_fail[0][0]}"]')
    lines.append("#[test]")
    if rel.as_posix() in expected:
        # (Passes while the body fails; fails once it passes: expected_failures.txt.)
        lines.append("#[should_panic]")
        lines.append(f"fn {fn_name}() {{")
        lines.append("    expected_failure(|| {")
        indent = "        "
    else:
        lines.append(f"fn {fn_name}() -> DrawResult<()> {{")
        indent = "    "
    lines.append(indent + "let mut d = Draw::new();")
    for stmt in res.rust_body:
        lines.append(indent + stmt)
    lines.append(indent + "Ok(())")
    if rel.as_posix() in expected:
        lines.append("    })")
    lines.append("}")
    return "\n".join(lines)


HEADER = "// @generated by tools/draw2rs.py from OCCT {src}. Do not edit; re-run the tool instead.\n"


def write_if_changed(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text() != text:
        path.write_text(text)


def existing_command_methods() -> set:
    names: set = set()
    for f in (TESTS_CRATE / "src" / "draw").rglob("*.rs"):
        names.update(re.findall(r"\bpub fn (r#)?([A-Za-z0-9_]+)\s*\(", f.read_text()))
    return {n[1] for n in names}


# Scripts not converted that the Tcl interpreter does not run either: data exchange, viewers,
# and those draw_deviations.json skips.
NOT_RUN_PREFIXES = ("data-exchange:", "viewer:", "draw_deviations.json:")


def vendor_scripts(occt: Path, skipped: list) -> set:
    """Copies the scripts not converted (that the Tcl interpreter can run) into `scripts/`,
    with their groups' and grids' begin and end files and parse.rules, for `tests/draw_tcl.rs`;
    the scripts' paths (`tests/...`), as converted."""
    import shutil
    root = TESTS_CRATE / "scripts"
    tests = occt / "tests"
    if root.exists():
        shutil.rmtree(root)
    root.mkdir()
    (root / "README.md").write_text(
        "# OCCT's DRAW scripts, as they are\n\n"
        "Copied from OCCT's `tests` by `tools/draw2rs.py`: the scripts it does not convert to Rust,\n"
        "which `tests/draw_tcl.rs` runs on a Tcl interpreter (`src/draw/tcl.rs`), with their groups'\n"
        "and grids' `begin` and `end` files and `parse.rules`. LGPL-2.1 with the Open CASCADE\n"
        "exception, as OCCT (see `../README.md`). Do not edit; re-run the tool instead.\n")
    copied = set()

    def copy(rel: Path) -> None:
        src = tests / rel
        if src.is_file() and rel not in copied:
            (root / rel).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, root / rel)
            copied.add(rel)

    copy(Path("parse.rules"))
    run = set()
    for path, reason in skipped:
        if reason.startswith(NOT_RUN_PREFIXES):
            continue
        rel = Path(path)
        group, grid = rel.parts[0], rel.parts[1]
        copy(rel)
        for d in (Path(group), Path(group) / grid):
            for f in ("begin", "end", "parse.rules"):
                copy(d / f)
        run.add(f"tests/{path}")
    return run


def write_dispatch() -> None:
    """`src/draw/dispatch.rs`: each DRAW command `Draw` has a method for, called by its name
    (for scripts run by the Tcl interpreter, see `src/draw/tcl.rs`)."""
    methods: set = set()
    for f in (TESTS_CRATE / "src" / "draw").rglob("*.rs"):
        methods.update(m[0] + m[1] for m in re.findall(r"\bpub fn (r#)?([A-Za-z0-9_]+)\(&mut self, args: &\[&str\]\) -> DrawResult", f.read_text()))
    arms = []
    for m in sorted(methods, key=lambda x: x.replace("r#", "")):
        name = m.replace("r#", "")
        if name.endswith("_cmd") and name[:-4] in RESERVED_METHODS:
            name = name[:-4]
        if re.fullmatch(r"_\d.*", name):
            name = name[1:]
        if cmd_method(name) != m:
            continue
        arms.append(f'        "{name}" => d.{m}(args),')
    text = (HEADER.format(src="commands (the methods of `Draw`)") + "\n"
            "use super::{Draw, DrawResult};\n\n"
            "/// Runs the DRAW command `name` with `args`; `None` if `Draw` has no such command.\n"
            "pub fn dispatch(d: &mut Draw, name: &str, args: &[&str]) -> Option<DrawResult> {\n"
            "    Some(match name {\n" + "\n".join(arms) + "\n        _ => return None,\n    })\n}\n")
    write_if_changed(TESTS_CRATE / "src" / "draw" / "dispatch.rs", text)


def occt_version(occt: Path) -> str:
    """OCCT's version, as `git describe` names its checkout (the report states no local path)."""
    try:
        out = subprocess.run(["git", "-C", str(occt), "describe", "--tags"], capture_output=True, text=True, check=True)
        return f"`{out.stdout.strip()}`"
    except (OSError, subprocess.CalledProcessError):
        return "unknown version"


def append_shims(commands: Counter, examples: dict, scripts_using: dict) -> list[str]:
    have = existing_command_methods()
    missing = sorted(c for c in commands if cmd_method(c).replace("r#", "") not in have)
    if not missing:
        return []
    out = ["", "// Shims appended by tools/draw2rs.py. Implement in place; the tool never",
           "// modifies an existing method.", "#[allow(non_snake_case)]", "impl Draw {"]
    for c in missing:
        m = cmd_method(c)
        ex = re.sub(r"\s+", " ", examples[c]).replace("`", "'")
        out.append(f"    /// DRAW command `{c}`: {commands[c]} uses in {scripts_using[c]} converted scripts.")
        out.append("    ///")
        out.append(f"    /// Example: `{ex}`")
        out.append(f"    pub fn {m}(&mut self, args: &[&str]) -> DrawResult {{")
        out.append(f"        self.not_implemented({rust_str(c)}, args)")
        out.append("    }")
    out.append("}")
    if not COMMANDS_RS.exists():
        COMMANDS_RS.parent.mkdir(parents=True, exist_ok=True)
        COMMANDS_RS.write_text("//! DRAW command implementations, one method per command.\n\n"
                               "use super::{Draw, DrawResult};\n")
    with COMMANDS_RS.open("a") as f:
        f.write("\n".join(out) + "\n")
    return missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--occt", type=Path, default=Path(os.environ.get("OCCT_ROOT", Path.home() / "OCCT")))
    ap.add_argument("--data-dir", type=Path, action="append", default=[],
                    help="extra directories searched for test data (e.g. the OCCT dataset)")
    ap.add_argument("--groups", default=",".join(DEFAULT_GROUPS))
    ap.add_argument("--no-shims", action="store_true", help="don't append command shims")
    args = ap.parse_args()

    occt: Path = args.occt.expanduser().resolve()
    groups = [g for g in args.groups.split(",") if g]
    data_dirs = [d.expanduser().resolve() for d in args.data_dir] + [Path.home() / "occt-data", occt / "data"]
    data_index = build_data_index(data_dirs)
    cases = discover(occt, groups)

    commands: Counter = Counter()
    examples: dict = {}
    scripts_using: Counter = Counter()
    skipped: Counter = Counter()
    skipped_all: list = []
    skipped_examples: dict = defaultdict(list)
    per_group: dict = defaultdict(lambda: defaultdict(list))
    stats = Counter()

    deviations = load_deviations()
    expected = load_expected_failures()
    converted: set = set()
    selectors = load_selectors()
    pre_cache: dict = {}
    for case in cases:
        key_g, key_grid = case.path.parent.parent, case.path.parent
        for key in (key_g, key_grid):
            if key not in pre_cache:
                pre_cache[key] = preamble_info(key / "begin")
        gconsts, gprocs, _, gdsets = pre_cache[key_g]
        cconsts, cprocs, ccomplex, cdsets = pre_cache[key_grid]
        consts = {**gconsts, **cconsts}
        procs = gprocs | cprocs
        stats["total"] += 1
        try:
            entry = deviations.get(case.path.relative_to(occt).as_posix())
            if isinstance(entry, dict):
                raise Unsupported(f"draw_deviations.json: {entry['skip']}")
            if ccomplex:
                raise Unsupported("grid begin file")
            res = translate_case(case, consts, procs, data_index, gdsets + cdsets)
            if not res.commands:
                raise Unsupported("no commands")
            if case.group in UNTRANSLATED_END and not any(c.startswith("check") for c in res.commands):
                raise Unsupported("checks only in the group's end script")
        except Unsupported as e:
            reason = str(e)
            skipped_all.append((case.path.relative_to(occt / "tests").as_posix(), reason))
            skipped[reason] += 1
            if len(skipped_examples[reason]) < 3:
                skipped_examples[reason].append(case.path.relative_to(occt / "tests").as_posix())
            continue
        stats["converted"] += 1
        if any(f not in data_index for f in res.data_files):
            stats["missing_data"] += 1
        commands.update(res.commands)
        for c in res.commands:
            scripts_using[c] += 1
            examples.setdefault(c, res.examples[c])
        per_group[case.group][case.grid].append((case, res))

    # write tests
    tests_dir = TESTS_CRATE / "tests"
    for old in tests_dir.glob("draw_*"):
        if old.is_dir() and old.name[len("draw_"):] not in per_group:
            for f in old.glob("*.rs"):
                f.unlink()
            old.rmdir()
    git_hash = subprocess.run(["git", "-C", str(occt), "rev-parse", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
    cases_dir = TESTS_CRATE / "cases" / "occt" / "draw"
    yaml_live: set = set()
    for group, grids in per_group.items():
        gname = re.sub(r"[^a-z0-9_]", "_", group.lower())
        mods = []
        gdir = tests_dir / f"draw_{gname}"
        live = set()
        for grid, items in sorted(grids.items()):
            mname = re.sub(r"[^A-Za-z0-9_]", "_", grid).lower()
            if mname[0].isdigit() or mname in RUST_KEYWORDS:
                mname = "g_" + mname
            used: set = set()
            body = [HEADER.format(src=f"tests/{group}/{grid}"), "use occt_tests::draw::*;", ""]
            docs = []
            for case, res in items:
                converted.add(case.path.relative_to(occt).as_posix())
                doc = try_yaml(case, res, data_index, occt, deviations, selectors, git_hash, expected)
                if doc is not None:
                    docs.append(doc)
                    stats["yaml"] += 1
                    continue
                body.append(emit_case(case, res, rust_test_name(case.name, used), data_index, occt,
                                      deviations, selectors, expected))
                body.append("")
            if docs:
                ypath = cases_dir / gname / f"{mname}.yaml"
                header = f"# Generated by tools/draw2rs.py from OCCT tests/{group}/{grid}. Do not edit; re-run the tool instead.\n"
                write_if_changed(ypath, header + "---\n".join(docs))
                yaml_live.add(ypath)
            if used:
                mods.append(mname)
                write_if_changed(gdir / f"{mname}.rs", "\n".join(body))
                live.add(f"{mname}.rs")
        for old in gdir.glob("*.rs") if gdir.exists() else []:
            if old.name not in live and old.name != "main.rs":
                old.unlink()
        if mods:
            root = [HEADER.format(src=f"tests/{group}"), "#![allow(unused_variables)]", ""]
            root += [f"mod {m};" for m in mods]
            write_if_changed(gdir / "main.rs", "\n".join(root) + "\n")
        elif gdir.exists():
            for f in gdir.glob("*.rs"):
                f.unlink()
            gdir.rmdir()
    for old in cases_dir.rglob("*.yaml") if cases_dir.exists() else []:
        if old not in yaml_live:
            old.unlink()

    added = [] if args.no_shims else append_shims(commands, examples, scripts_using)
    write_dispatch()
    # (The scripts not converted, with why: for running them on the Tcl interpreter.)
    target = TESTS_CRATE.parent.parent / "target"
    target.mkdir(exist_ok=True)
    (target / "draw_skipped.txt").write_text("".join(f"{p}\t{r}\n" for p, r in skipped_all))
    converted |= vendor_scripts(occt, skipped_all)

    stale = sorted(expected - converted)
    if stale:
        sys.exit(f"{EXPECTED_FAILURES.name}: not converted any more (remove them): {', '.join(stale)}")

    # report
    rep = [f"# DRAW conversion report", "",
           f"OCCT: {occt_version(occt)}; groups: {', '.join(groups)}", "",
           f"- scripts considered: {stats['total']}",
           f"- converted: {stats['converted']} ({stats['missing_data']} ignored for missing data), "
           f"{stats['yaml']} of them as YAML cases (cases/occt/draw)",
           f"- skipped: {stats['total'] - stats['converted']}",
           f"- distinct DRAW commands used by converted scripts: {len(commands)}", "",
           "## Skip reasons", "", "| reason | scripts | examples |", "|---|---|---|"]
    for r, n in skipped.most_common():
        rep.append(f"| {r} | {n} | {', '.join(skipped_examples[r])} |")
    rep += ["", "## Commands by use", "", "| command | uses | scripts |", "|---|---|---|"]
    for c, n in commands.most_common():
        rep.append(f"| `{c}` | {n} | {scripts_using[c]} |")
    write_if_changed(TESTS_CRATE / "DRAW_REPORT.md", "\n".join(rep) + "\n")

    print(f"converted {stats['converted']}/{stats['total']} scripts ({stats['yaml']} as YAML) "
          f"({stats['missing_data']} need missing data); {len(commands)} commands; "
          f"{len(added)} new shims")
    for r, n in skipped.most_common(12):
        print(f"  skipped {n:5d}  {r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
