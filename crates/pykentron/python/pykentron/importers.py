"""Reading shapes from files (CadQuery's `importers`)."""

from __future__ import annotations

from . import _kentron as native
from .geom import Color, Location
from .shapes import Shape


def importBrep(fileName) -> "Workplane":
    from .workplane import Workplane

    return Workplane("XY").add(Shape.importBrep(str(fileName)))


def importStep(fileName, unit: str = "MM") -> "Workplane":
    """The file's shapes (one compound of its products, placed) on a workplane."""
    from .workplane import Workplane

    if unit.upper() != "MM":
        raise NotImplementedError("STEP units other than millimetres")
    return Workplane("XY").add(Shape.cast(native.read_step(str(fileName))))


def importStepAssembly(fileName):
    """The file's product structure as an assembly: occurrences named and placed, parts with
    their shapes and colours."""
    from .assembly import Assembly

    def build(node):
        name, handle, rows, color, children = node
        a = Assembly(
            Shape.cast(handle) if handle is not None else None,
            loc=Location([list(r) for r in rows]),
            name=name,
            color=Color(*color) if color is not None else None,
        )
        for child in children:
            a.add(build(child))
        return a

    roots = native.read_step_assembly(str(fileName))
    if len(roots) == 1:
        return build(roots[0])
    top = Assembly(name="root")
    for r in roots:
        top.add(build(r))
    return top


def importDXF(*args, **kwargs):
    raise NotImplementedError("reading DXF files")
