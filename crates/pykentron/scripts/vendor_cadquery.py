#!/usr/bin/env python3
"""Copies CadQuery's pure-Python layers (no OCCT calls) into the `cadquery` shim, unchanged but
for the few patches in PATCHES (imports deferred to first use, for a quicker import).

CadQuery's `Workplane`, `Sketch`, selectors, `cqgi`, free functions and assemblies are written
on its `occ_impl` layer, which the shim provides from pykentron; this copies those modules from
a CadQuery release (Apache-2.0: its licence is copied with them). Run it again to update them
after changing the release.

Usage:
    python3 crates/pykentron/scripts/vendor_cadquery.py [--source DIR]
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHIM = HERE.parent / "python" / "pykentron" / "cadquery_shim" / "cadquery"
EXPORTERS = HERE.parent / "python" / "pykentron" / "occ_impl" / "exporters"
REPOSITORY = "https://github.com/CadQuery/cadquery.git"
RELEASE = "v2.8.0"

# The modules copied (relative to CadQuery's `cadquery` package): those with no direct OCCT
# calls; `__init__.py` and the rest of `occ_impl` are the shim's own.
VENDORED = [
    "assembly.py",
    "cq.py",
    "cqgi.py",
    "func.py",
    "hull.py",
    "selectors.py",
    "sketch.py",
    "types.py",
    "units.py",
    "utils.py",
    "plugins/__init__.py",
    "occ_impl/shape_protocols.py",
    "occ_impl/sketch_solver.py",
    "occ_impl/solver.py",
]

# Changes to the copies, `(module, CadQuery's text, replacement)`: what the import would
# otherwise load or build up front, deferred to first use (numpy for hulls and the sketch
# solver, nlopt, pyparsing and the selector grammars for string selectors: half the shim's
# import). Each text must occur
# once, or vendoring stops.
PATCHES = [
    (
        "hull.py",
        "from numpy import lexsort, argmin, argmax\n",
        "# (pykentron: numpy imported on first use, not with the module.)\n"
        "def lexsort(*args, **kwargs):\n"
        "    from numpy import lexsort\n\n"
        "    return lexsort(*args, **kwargs)\n\n\n"
        "def argmin(*args, **kwargs):\n"
        "    from numpy import argmin\n\n"
        "    return argmin(*args, **kwargs)\n\n\n"
        "def argmax(*args, **kwargs):\n"
        "    from numpy import argmax\n\n"
        "    return argmax(*args, **kwargs)\n\n",
    ),
    (
        "occ_impl/sketch_solver.py",
        "from numpy.typing import NDArray as Array\n"
        "from numpy import float64 as Float\n",
        "# (pykentron: numpy and nlopt imported on first use, not with the module: the names\n"
        "# below stand in for theirs.)\n"
        "class Array:\n"
        "    \"\"\"numpy.typing.NDArray, in annotations.\"\"\"\n\n"
        "    def __class_getitem__(cls, item):\n"
        "        return cls\n\n\n"
        "Float = float\n",
    ),
    (
        "occ_impl/sketch_solver.py",
        "from numpy import array, full, inf, sign\n"
        "from numpy.linalg import norm\n"
        "import nlopt\n",
        "inf = float(\"inf\")\n\n\n"
        "def array(*args, **kwargs):\n"
        "    from numpy import array\n\n"
        "    return array(*args, **kwargs)\n\n\n"
        "def full(*args, **kwargs):\n"
        "    from numpy import full\n\n"
        "    return full(*args, **kwargs)\n\n\n"
        "def sign(*args, **kwargs):\n"
        "    from numpy import sign\n\n"
        "    return sign(*args, **kwargs)\n\n\n"
        "def norm(*args, **kwargs):\n"
        "    from numpy.linalg import norm\n\n"
        "    return norm(*args, **kwargs)\n\n\n"
        "class _OnFirstUse:\n"
        "    \"\"\"A module imported when one of its names is first read.\"\"\"\n\n"
        "    def __init__(self, name):\n"
        "        self._name = name\n\n"
        "    def __getattr__(self, attr):\n"
        "        import importlib\n\n"
        "        return getattr(importlib.import_module(self._name), attr)\n\n\n"
        "nlopt = _OnFirstUse(\"nlopt\")\n",
    ),
    (
        "selectors.py",
        "from pyparsing import (\n"
        "    pyparsing_common,\n"
        "    Literal,\n"
        "    Word,\n"
        "    nums,\n"
        "    Optional,\n"
        "    Combine,\n"
        "    one_of,\n"
        "    Group,\n"
        "    infix_notation,\n"
        "    opAssoc,\n"
        ")\n",
        "# (pykentron: pyparsing imported, and the grammars below made, on first use.)\n",
    ),
    (
        "selectors.py",
        "_grammar = _makeGrammar()  # make a grammar instance\n",
        "",
    ),
    (
        "selectors.py",
        "_expression_grammar = _makeExpressionGrammar(_grammar)\n",
        "def _grammars():\n"
        "    \"\"\"(pykentron) The selector grammars, made on first use.\"\"\"\n"
        "    if \"_expression_grammar\" not in globals():\n"
        "        import pyparsing\n\n"
        "        for name in (\"pyparsing_common\", \"Literal\", \"Word\", \"nums\", \"Optional\", \"Combine\", \"one_of\", \"Group\", \"infix_notation\", \"opAssoc\"):\n"
        "            globals()[name] = getattr(pyparsing, name)\n"
        "        globals()[\"_grammar\"] = _makeGrammar()\n"
        "        globals()[\"_expression_grammar\"] = _makeExpressionGrammar(_grammar)\n"
        "    return _grammar, _expression_grammar\n\n\n"
        "def __getattr__(name):\n"
        "    if name in (\"_grammar\", \"_expression_grammar\"):\n"
        "        return _grammars()[name == \"_expression_grammar\"]\n"
        "    raise AttributeError(f\"module {__name__!r} has no attribute {name!r}\")\n",
    ),
    (
        "selectors.py",
        "        parse_result = _expression_grammar.parse_string(selectorString, parse_all=True)\n",
        "        parse_result = _grammars()[1].parse_string(selectorString, parse_all=True)\n",
    ),
]


def patch(path: Path, old: str, new: str) -> None:
    text = path.read_text()
    if text.count(old) != 1:
        raise SystemExit(f"{path.name}: the text to patch occurs {text.count(old)} times, not once: {old[:60]!r}")
    path.write_text(text.replace(old, new))


# CadQuery's mesh writers (pure Python, on shapes' tessellations), into pykentron's
# `occ_impl.exporters`, where its other exporters are.
VENDORED_EXPORTERS = ["amf.py", "json.py", "threemf.py"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", help="a CadQuery checkout at the release (fetched if not given)")
    args = ap.parse_args()
    if args.source:
        src = Path(args.source)
    else:
        src = HERE.parent.parent.parent / "target" / "cadquery-tests" / f"cadquery-{RELEASE}"
        if not (src / "cadquery").is_dir():
            subprocess.run(["git", "clone", "--quiet", "--depth", "1", "--branch", RELEASE, REPOSITORY, str(src)], check=True)
    for rel in VENDORED:
        target = SHIM / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src / "cadquery" / rel, target)
    for rel, old, new in PATCHES:
        patch(SHIM / rel, old, new)
    for rel in VENDORED_EXPORTERS:
        shutil.copyfile(src / "cadquery" / "occ_impl" / "exporters" / rel, EXPORTERS / rel)
    for target, files, prefix in ((SHIM, VENDORED, ""), (EXPORTERS, VENDORED_EXPORTERS, "occ_impl/exporters/")):
        shutil.copyfile(src / "LICENSE", target / "LICENSE.CadQuery")
        patched = sorted({rel for rel, _, _ in PATCHES}) if target == SHIM else []
        (target / "VENDORED.txt").write_text(
            f"Copied unchanged from CadQuery {RELEASE} ({REPOSITORY}), Apache License 2.0 "
            "(LICENSE.CadQuery), by crates/pykentron/scripts/vendor_cadquery.py:\n" + "".join(f"  {prefix}{r}\n" for r in files)
            + (f"Patched there (imports deferred to first use): {', '.join(patched)}.\n" if patched else "")
        )
    print(f"copied {len(VENDORED) + len(VENDORED_EXPORTERS)} modules from CadQuery {RELEASE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
