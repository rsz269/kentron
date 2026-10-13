"""CadQuery's `occ_impl.importers` on pykentron: files read onto CadQuery's workplanes."""

from __future__ import annotations

import sys
from typing import List

from ... import _kentron as native
from ...shapes import Shape


class ImportTypes:
    STEP = "STEP"
    DXF = "DXF"
    BREP = "BREP"
    BIN = "BIN"


class UNITS:
    MM = "mm"
    IN = "in"


# Millimetres per unit (CadQuery's unit names); kentron reads STEP files in millimetres.
_MILLIMETRES = {"MM": 1.0, "CM": 10.0, "M": 1000.0, "KM": 1e6, "INCH": 25.4, "FT": 304.8, "MI": 1609344.0, "UM": 1e-3, "NM": 1e-6}


def _workplane():
    """CadQuery's Workplane where the `cadquery` shim runs it, else pykentron's."""
    cq = sys.modules.get("cadquery.cq")
    if cq is not None:
        return cq.Workplane
    from ...workplane import Workplane

    return Workplane


def importShape(importType: str, fileName: str, unit: str = "MM", *args, **kwargs):
    if importType == ImportTypes.STEP:
        return importStep(fileName, unit)
    if importType == ImportTypes.DXF:
        return importDXF(fileName, *args, **kwargs)
    if importType == ImportTypes.BREP:
        return importBrep(fileName)
    if importType == ImportTypes.BIN:
        return importBin(fileName)
    raise RuntimeError(f"Unsupported import type: {importType!r}")


def importBrep(fileName: str):
    """The BREP file's shape (one, a compound if several) on a workplane."""
    return _workplane()("XY").newObject([Shape.importBrep(str(fileName))])


def importBin(fileName: str):
    raise NotImplementedError("reading OCCT's binary BREP files")


def importStep(fileName: str, unit: str = "MM"):
    """The STEP file's shapes on a workplane, in `unit`s."""
    try:
        factor = _MILLIMETRES[unit.upper()]
    except KeyError:
        raise ValueError(f"unknown unit {unit!r}")
    from ...importers import _step_roots

    try:
        shapes = _step_roots(fileName)
    except (OSError, ValueError) as e:
        raise ValueError("STEP File could not be loaded") from e
    # A shape per root product (as OCCT's reader), in `unit`s.
    if factor != 1.0:
        shapes = [s.scale(1.0 / factor) for s in shapes]
    return _workplane()("XY").newObject(shapes)


def importDXF(filename: str, tol: float = 1e-6, exclude: List[str] = [], include: List[str] = []):
    """The faces of a DXF file's (selected) layers on a workplane."""
    from .dxf import _importDXF

    return _workplane()("XY").newObject(_importDXF(filename, tol, exclude, include))
