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
    "occ_impl/nurbs.py",
    "occ_impl/jupyter_tools.py",
    "vis.py",
]
# (Copied under another name: CadQuery's rendering, which needs VTK, behind the shim's own
# `vis`, which takes it where VTK is installed.)
VENDORED_AS = {"vis.py": "_vis_vtk.py"}

# Changes to the copies, `(module, CadQuery's text, replacement)`: what the import would
# otherwise load or build up front, deferred to first use (numpy for hulls and the sketch
# solver, nlopt, casadi, pyparsing and the selector grammars, which string selectors no longer
# need: they are read by pykentron.selector_syntax), half the shim's import; and the two
# constraint solvers' classes handed to pykentron.constraint_solvers, which solves their
# problems with kentron (no casadi, nlopt or numpy) unless PYKENTRON_SOLVER=cadquery. Each
# text must occur once, or vendoring stops.
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
        "occ_impl/solver.py",
        "from typing import (\n",
        "from __future__ import annotations  # (pykentron: annotations not evaluated, so casadi's\n"
        "# names in them do not import it.)\n\n"
        "from typing import (\n",
    ),
    (
        "occ_impl/solver.py",
        "import casadi as ca\n",
        "class _OnFirstUse:\n"
        "    \"\"\"(pykentron) A module imported when one of its names is first read.\"\"\"\n\n"
        "    def __init__(self, name):\n"
        "        self._name = name\n\n"
        "    def __getattr__(self, attr):\n"
        "        import importlib\n\n"
        "        return getattr(importlib.import_module(self._name), attr)\n\n\n"
        "ca = _OnFirstUse(\"casadi\")\n",
    ),
    (
        "occ_impl/solver.py",
        "        return locs, result\n",
        "        return locs, result\n\n\n"
        "# (pykentron: the problem solved by kentron, unless PYKENTRON_SOLVER=cadquery.)\n"
        "from pykentron.constraint_solvers import assembly_solver as _assembly_solver\n\n"
        "ConstraintSolver = _assembly_solver(ConstraintSolver)\n",
    ),
    (
        "occ_impl/sketch_solver.py",
        "        return [x[i1:i2] for i1, i2 in zip(ixs, ixs[1:])], status\n",
        "        return [x[i1:i2] for i1, i2 in zip(ixs, ixs[1:])], status\n\n\n"
        "# (pykentron: the problem solved by kentron, unless PYKENTRON_SOLVER=cadquery; with it,\n"
        "# arcs' points as tuples, without numpy.)\n"
        "from pykentron.constraint_solvers import sketch_solver as _sketch_solver\n\n"
        "SketchConstraintSolver, arc_first, arc_last, arc_point = _sketch_solver(SketchConstraintSolver, arc_first, arc_last, arc_point)\n",
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
        "def _parse(text):\n"
        "    \"\"\"(pykentron) The selector `text` describes, read as the grammars above read it, by\n"
        "    pykentron.selector_syntax (without pyparsing, which only the grammars' own users load).\"\"\"\n"
        "    from pykentron.selector_syntax import parse\n\n"
        "    return parse(\n"
        "        text,\n"
        "        set(geom_LUT_EDGE.values()) | set(geom_LUT_FACE.values()),\n"
        "        _SimpleStringSyntaxSelector,\n"
        "        lambda items: reduce(AndSelector, items),\n"
        "        lambda items: reduce(SumSelector, items),\n"
        "        lambda items: reduce(SubtractSelector, items),\n"
        "        InverseSelector,\n"
        "    )\n\n\n"
        "def __getattr__(name):\n"
        "    if name in (\"_grammar\", \"_expression_grammar\"):\n"
        "        return _grammars()[name == \"_expression_grammar\"]\n"
        "    raise AttributeError(f\"module {__name__!r} has no attribute {name!r}\")\n",
    ),
    (
        "selectors.py",
        "        parse_result = _expression_grammar.parse_string(selectorString, parse_all=True)\n"
        "        self.mySelector = parse_result.asList()[0]\n",
        "        self.mySelector = _parse(selectorString)\n",
    ),
]


# CadQuery's NURBS module (numpy, scipy and numba arithmetic on B-splines) and its rendering
# (VTK), which reach OCCT only to turn their B-splines into edges and faces and back: those
# conversions kentron's, through pykentron's native module.
PATCHES += [
    (
        "occ_impl/nurbs.py",
        "from OCP.Geom import Geom_BSplineCurve, Geom_BSplineSurface\n"
        "from OCP.TColgp import TColgp_Array1OfPnt, TColgp_Array2OfPnt\n"
        "from OCP.TColStd import (\n"
        "    TColStd_Array1OfInteger,\n"
        "    TColStd_Array1OfReal,\n"
        ")\n"
        "from OCP.gp import gp_Pnt\n"
        "from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace\n",
        "# (pykentron: no OCCT. Its B-splines the OCP shim's, as data; the conversions to\n"
        "# edges and faces and back kentron's.)\n"
        "from OCP.Geom import Geom_BSplineCurve, Geom_BSplineSurface\n"
        "TColgp_Array1OfPnt = TColgp_Array2OfPnt = TColStd_Array1OfInteger = TColStd_Array1OfReal = object\n"
        "from pykentron import _kentron as _native\n",
    ),
    (
        "occ_impl/nurbs.py",
        "    def curve(self) -> Geom_BSplineCurve:\n\n"
        "        if self.periodic:\n",
        "    def _knots(self):\n"
        "        # (pykentron: as OCCT's knots and multiplicities, below.)\n"
        "        if self.periodic:\n"
        "            return list(self.knots), [1] * len(self.knots)\n"
        "        unique_knots, mults_arr = np.unique(self.knots, return_counts=True)\n"
        "        return list(unique_knots), [int(m) for m in mults_arr]\n\n"
        "    def curve(self) -> Geom_BSplineCurve:\n\n"
        "        knots, mults = self._knots()\n"
        "        return Geom_BSplineCurve(self.pts, knots, mults, self.order, self.periodic)\n"
        "        if self.periodic:\n",
    ),
    (
        "occ_impl/nurbs.py",
        "        return Edge(BRepBuilderAPI_MakeEdge(self.curve()).Shape())\n",
        "        knots, mults = self._knots()\n"
        "        poles = [tuple(map(float, p)) for p in self.pts]\n"
        "        return Edge(_native.make_bspline(poles, None, knots, mults, self.order, self.periodic))\n",
    ),
    (
        "occ_impl/nurbs.py",
        "        g = e._geomAdaptor().BSpline()\n\n"
        "        knots = np.repeat(list(g.Knots()), list(g.Multiplicities()))\n"
        "        pts = np.array([(p.X(), p.Y(), p.Z()) for p in g.Poles()])\n"
        "        order = g.Degree()\n"
        "        periodic = g.IsPeriodic()\n",
        "        order, poles, knots_, mults, _, periodic = _native.edge_bspline_data(e.wrapped)\n\n"
        "        knots = np.repeat(knots_, mults)\n"
        "        pts = np.array(poles)\n",
    ),
    (
        "occ_impl/nurbs.py",
        "        Convert to a OCCT Surface object.\n"
        "        \"\"\"\n",
        "        Convert to a OCCT Surface object.\n"
        "        \"\"\"\n\n"
        "        uknots, umults = np.unique(self.uknots, return_counts=True)\n"
        "        vknots, vmults = np.unique(self.vknots, return_counts=True)\n"
        "        return Geom_BSplineSurface(self.pts, uknots, vknots, umults, vmults, self.uorder, self.vorder, self.uperiodic, self.vperiodic)\n",
    ),
    (
        "occ_impl/nurbs.py",
        "        return Face(BRepBuilderAPI_MakeFace(self.surface(), tol).Shape())\n",
        "        uknots, umults = np.unique(self.uknots, return_counts=True)\n"
        "        vknots, vmults = np.unique(self.vknots, return_counts=True)\n"
        "        poles = [[tuple(map(float, p)) for p in row] for row in self.pts]\n"
        "        return Face(\n"
        "            _native.make_bspline_face(\n"
        "                poles, None, list(uknots), [int(m) for m in umults], list(vknots), [int(m) for m in vmults],\n"
        "                self.uorder, self.vorder, self.uperiodic, self.vperiodic,\n"
        "            )\n"
        "        )\n",
    ),
    (
        "occ_impl/nurbs.py",
        "        g = cast(Geom_BSplineSurface, f._geomAdaptor())\n\n"
        "        uknots = np.repeat(list(g.UKnots()), list(g.UMultiplicities()))\n"
        "        vknots = np.repeat(list(g.VKnots()), list(g.VMultiplicities()))\n\n"
        "        tmp = []\n"
        "        for i in range(1, g.NbUPoles() + 1):\n"
        "            tmp.append(\n"
        "                [\n"
        "                    [g.Pole(i, j).X(), g.Pole(i, j).Y(), g.Pole(i, j).Z(),]\n"
        "                    for j in range(1, g.NbVPoles() + 1)\n"
        "                ]\n"
        "            )\n\n"
        "        pts = np.array(tmp)\n\n"
        "        uorder = g.UDegree()\n"
        "        vorder = g.VDegree()\n\n"
        "        uperiodic = g.IsUPeriodic()\n"
        "        vperiodic = g.IsVPeriodic()\n",
        "        (uorder, vorder), poles, (uk, um), (vk, vm), _, (uperiodic, vperiodic) = _native.face_bspline_data(f.wrapped)\n\n"
        "        uknots = np.repeat(uk, um)\n"
        "        vknots = np.repeat(vk, vm)\n"
        "        pts = np.array(poles)\n",
    ),
    (
        "vis.py",
        "from OCP.TopoDS import TopoDS_Shape\n"
        "from OCP.Geom import Geom_BSplineSurface\n",
        "# (pykentron: no OCCT shapes besides its own; ctrlPts reads B-splines through nurbs.)\n"
        "class TopoDS_Shape:\n"
        "    pass\n",
    ),
    (
        "vis.py",
        "        if isinstance(s._geomAdaptor(), Geom_BSplineSurface):\n"
        "            surf = cast(Geom_BSplineSurface, s._geomAdaptor())\n"
        "        else:\n"
        "            raise ValueError(\n"
        "                f\"Only NURBS surfaces are supported, encountered {s._geomAdaptor()}\"\n"
        "            )\n\n"
        "        Nu = surf.NbUPoles()\n"
        "        Nv = surf.NbVPoles()\n\n"
        "        u_periodic = surf.IsUPeriodic()\n"
        "        v_periodic = surf.IsVPeriodic()\n\n"
        "        # add points\n"
        "        for i in range(Nu):\n"
        "            for j in range(Nv):\n"
        "                pt = surf.Pole(i + 1, j + 1)\n"
        "                points.InsertNextPoint(pt.X(), pt.Y(), pt.Z())\n",
        "        if s.geomType() != \"BSPLINE\":\n"
        "            raise ValueError(f\"Only NURBS surfaces are supported, encountered {s.geomType()}\")\n"
        "        surf = Surface.fromFace(s)\n\n"
        "        Nu, Nv = surf.pts.shape[:2]\n"
        "        u_periodic, v_periodic = surf.uperiodic, surf.vperiodic\n\n"
        "        # add points\n"
        "        for i in range(Nu):\n"
        "            for j in range(Nv):\n"
        "                points.InsertNextPoint(*map(float, surf.pts[i, j]))\n",
    ),
    (
        "vis.py",
        "        if s.geomType() == \"BSPLINE\":\n"
        "            curve = s._geomAdaptor().BSpline()\n",
        "        if s.geomType() == \"BSPLINE\":\n"
        "            curve = Curve.fromEdge(s)\n",
    ),
    (
        "vis.py",
        "        for pt in curve.Poles():\n"
        "            points.InsertNextPoint(pt.X(), pt.Y(), pt.Z())\n\n"
        "        N = curve.NbPoles()\n",
        "        for pt in curve.pts:\n"
        "            points.InsertNextPoint(*map(float, pt))\n\n"
        "        N = len(curve.pts)\n",
    ),
    (
        "vis.py",
        "        if curve.IsPeriodic():\n",
        "        if curve.periodic:\n",
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
        target = SHIM / VENDORED_AS.get(rel, rel)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src / "cadquery" / rel, target)
    for rel, old, new in PATCHES:
        patch(SHIM / VENDORED_AS.get(rel, rel), old, new)
    for rel in VENDORED_EXPORTERS:
        shutil.copyfile(src / "cadquery" / "occ_impl" / "exporters" / rel, EXPORTERS / rel)
    for target, files, prefix in ((SHIM, VENDORED, ""), (EXPORTERS, VENDORED_EXPORTERS, "occ_impl/exporters/")):
        shutil.copyfile(src / "LICENSE", target / "LICENSE.CadQuery")
        patched = sorted({rel for rel, _, _ in PATCHES}) if target == SHIM else []
        (target / "VENDORED.txt").write_text(
            f"Copied unchanged from CadQuery {RELEASE} ({REPOSITORY}), Apache License 2.0 "
            "(LICENSE.CadQuery), by crates/pykentron/scripts/vendor_cadquery.py:\n" + "".join(f"  {prefix}{r}" + (f" (as {VENDORED_AS[r]})" if r in VENDORED_AS and target == SHIM else "") + "\n" for r in files)
            + (f"Patched there (imports deferred to first use; the constraint solvers kentron's; B-splines made and read by kentron, not OCCT): {', '.join(patched)}.\n" if patched else "")
        )
    print(f"copied {len(VENDORED) + len(VENDORED_EXPORTERS)} modules from CadQuery {RELEASE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
