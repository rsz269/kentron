"""Reading shapes from files (CadQuery's `importers`)."""

from __future__ import annotations

from . import _kentron as native
from .geom import Color, Location
from .shapes import Compound, Shape


def importBrep(fileName) -> "Workplane":
    from .workplane import Workplane

    return Workplane("XY").add(Shape.importBrep(str(fileName)))


def _step_roots(fileName) -> list:
    """The STEP file's shapes, one per root product (as OCCT's reader transfers its roots),
    in millimetres: a product's own shape and its child products, placed (a compound of
    them if more than one)."""
    from .geom import Location

    def shape_of(node):
        name, handle, rows, rgb, children = node
        parts = [Shape.cast(handle)] if handle is not None else []
        for child in children:
            part = shape_of(child)
            if part is not None:
                parts.append(part.moved(Location([list(r) for r in child[2]])))
        if not parts:
            return None
        return parts[0] if len(parts) == 1 else Compound.makeCompound(parts)

    roots = []
    for root in native.read_step_assembly(str(fileName)):
        shape = shape_of(root)
        if shape is not None:
            roots.append(shape.moved(Location([list(r) for r in root[2]])))
    return roots


def importStep(fileName, unit: str = "MM") -> "Workplane":
    """The file's shapes (one compound of its products, placed) on a workplane."""
    from .workplane import Workplane

    from .exporters import _mm_per

    # Read in millimetres, then in `unit`: a shape per root product.
    factor = 1.0 / _mm_per(unit)
    shapes = [s.scale(factor) if factor != 1.0 else s for s in _step_roots(fileName)]
    return Workplane("XY").add(shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes))


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
