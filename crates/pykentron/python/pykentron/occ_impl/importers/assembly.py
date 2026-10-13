"""CadQuery's `occ_impl.importers.assembly` on pykentron: assemblies read from STEP files."""

from __future__ import annotations

from ... import _kentron as native
from ...geom import Color, Location
from ...shapes import Shape


def importStep(assy, path: str, unit: str = "MM") -> None:
    """Fills `assy` with the product structure of a STEP file (as CadQuery's): the top-level
    assembly's name and placement on `assy`, sub-assemblies and parts below it, named, placed
    and coloured."""
    from ...exporters import _mm_per

    roots = native.read_step_assembly(str(path))
    if not roots:
        raise ValueError(f"Error reading STEP file: {path}")
    # Read in millimetres: in `unit`, shapes and placements' translations scaled.
    factor = 1.0 / _mm_per(unit)
    if factor != 1.0:

        def scaled(node):
            name, handle, rows, rgb, children = node
            if handle is not None:
                handle = Shape.cast(handle).scale(factor).wrapped
            rows = [list(r[:3]) + [r[3] * factor] for r in rows]
            return (name, handle, rows, rgb, [scaled(c) for c in children])

        roots = [scaled(r) for r in roots]

    def leaf(node, parent) -> None:
        name, handle, rows, rgb, children = node
        shape = Shape.cast(handle)
        loc = Location([list(r) for r in rows])
        color = Color(*rgb) if rgb is not None else None
        # (CadQuery writes an assembly's own shape as a "<name>_part" child.)
        if name.endswith("_part"):
            parent.obj, parent.loc, parent.color = shape, loc, color
        else:
            parent.add(assy.__class__(shape, loc=loc, name=name, color=color))

    def fill(node, parent) -> None:
        for child in node[4]:
            name, handle, rows, rgb, children = child
            if children:
                sub = assy.__class__(name=name)
                fill(child, sub)
                parent.add(sub, loc=Location([list(r) for r in rows]), name=name, color=Color(*rgb) if rgb is not None else None)
            elif handle is not None:
                leaf(child, parent)

    top = roots[0]
    name, handle, rows, rgb, children = top
    if children:
        # The top-level assembly: its name and placement taken by `assy`.
        assy.objects.pop(assy.name)
        assy.name = name
        assy.objects[assy.name] = assy
        if rgb is not None:
            assy.color = Color(*rgb)
        assy.loc = Location([list(r) for r in rows])
        fill(top, assy)
    else:
        for root in roots:
            if root[1] is not None:
                leaf(root, assy)


def importXbf(assy, path: str) -> None:
    raise NotImplementedError("reading XBF files (OCAF's binary format)")


def importXml(assy, path: str) -> None:
    raise NotImplementedError("reading XML files (OCAF's XML format)")
