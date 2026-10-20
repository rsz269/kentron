"""CadQuery's `occ_impl.exporters.assembly` on pykentron: assemblies written to STEP files."""

from __future__ import annotations

from typing import Literal

from ... import _kentron as native
from ...shapes import Compound, Shape

class ExportModes:
    DEFAULT = "default"
    FUSED = "fused"


STEPExportModeLiterals = Literal["default", "fused"]


def _shapes_of(obj) -> list:
    if obj is None:
        return []
    if isinstance(obj, Shape):
        return [obj]
    if hasattr(obj, "vals"):
        return [o for o in obj.vals() if isinstance(o, Shape)]
    return []


def _subshapes(assy) -> list:
    """The node's faces, edges and wires with names, colours or layers of their own, as the
    native writer takes them: `(shape, name, (r, g, b, a), layer)`."""
    marks = {}
    for attribute, slot in (("_subshape_names", 1), ("_subshape_colors", 2), ("_subshape_layers", 3)):
        for s, value in getattr(assy, attribute, {}).items():
            entry = marks.setdefault(id(s), [s, None, None, None])
            entry[slot] = tuple(value.toTuple()) if slot == 2 else value
    return [(s.wrapped, name, color, layer) for s, name, color, layer in marks.values()]


def _material(assy):
    m = getattr(assy, "material", None)
    return None if m is None else (m.name, m.description, float(m.density), m.densityUnit)


def _tree(assy, factor: float = 1.0):
    """The assembly as the native writer's nested `(name, shape, placement rows, rgba,
    children, extras)` tuples, lengths scaled by `factor` (shapes and placements'
    translations); the extras its subshapes' names, colours and layers and its material."""
    shapes = _shapes_of(assy.obj)
    shape = None
    if shapes:
        whole = shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)
        shape = (whole.scale(factor) if factor != 1.0 else whole).wrapped
    rgba = tuple(assy.color.toTuple()) if assy.color is not None else None
    rows = [list(r[:3]) + [r[3] * factor] for r in assy.loc.rows]
    # (Subshapes are found by identity: on a scaled copy, none.)
    extras = {"subshapes": _subshapes(assy) if factor == 1.0 else [], "material": _material(assy)}
    return (assy.name, shape, rows, rgba, [_tree(ch, factor) for ch in assy.children], extras)


def _fused_tree(assy, glue: bool = False, tol=None):
    """The assembly's parts, placed, fused into one shape (as CadQuery's `toFusedCAF`), as one
    product of the assembly's name: each face of the fused shape coloured as the part it
    lies on (the part with a face nearest it)."""
    from ...shapes import Vertex

    parts = [(shape.moved(loc), color) for shape, _, loc, color in assy]
    if not parts:
        raise ValueError(f"Error: Assembly {assy.name} has no shapes.")
    if len(parts) == 1:
        fused = parts[0][0]
        if fused.ShapeType() == "Compound":
            fused = fused.fuse(glue=glue, tol=tol) if len(fused.Solids()) > 1 else fused
    else:
        fused = parts[0][0].fuse(*[p for p, _ in parts[1:]], glue=glue, tol=tol)
    subshapes = []
    coloured = [(p, c) for p, c in parts if c is not None]
    if coloured:
        for face in fused.Faces():
            at = face.Center()
            point = Vertex.makeVertex(at.x, at.y, at.z)
            _, color = min(coloured, key=lambda pc: min(f.distance(point) for f in pc[0].Faces()))
            subshapes.append((face.wrapped, None, tuple(color.toTuple()), None))
    identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
    return (assy.name, fused.wrapped, identity, None, [], {"subshapes": subshapes, "material": None})


def exportAssembly(assy, path: str, mode: str = "default", unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    """Writes the assembly to a STEP file, its tree as products and occurrences; fused, as
    one product of all its parts fused, its faces coloured as theirs."""
    import os

    from ...exporters import _mm_per

    # (As OCCT's writer: nothing written to a directory that does not exist, and False.)
    folder = os.path.dirname(os.path.abspath(str(path)))
    if not os.path.isdir(folder):
        return False
    out = (outputUnit or unit).upper()
    if mode == "fused":
        if _mm_per(unit) != _mm_per(out):
            raise NotImplementedError("fused assemblies written in another unit")
        native.write_step_assembly(str(path), [_fused_tree(assy, glue=bool(kwargs.get("glue", False)), tol=kwargs.get("fuzzy_tol"))], unit=out)
        return True
    tree = _tree(assy, _mm_per(unit) / _mm_per(out))
    # (STEP places no root: a top level placed otherwise than at the origin is written as
    # the one child, so placed, of a root of its name, which reading takes back off.)
    if any(abs(v - (1.0 if i == j else 0.0)) > 1e-12 for i, row in enumerate(assy.loc.rows) for j, v in enumerate(row)):
        identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        tree = (assy.name, None, identity, None, [tree], {"subshapes": [], "material": None})
    native.write_step_assembly(str(path), [tree], unit=out)
    return True


def exportStepMeta(assy, path: str, unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    return exportAssembly(assy, path, unit=unit, outputUnit=outputUnit)


def exportCAF(assy, path: str, binary: bool = True) -> bool:
    raise NotImplementedError(f"exporting assemblies to {'XBF' if binary else 'XML'} files")


def exportVRML(assy, path: str, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> None:
    raise NotImplementedError("exporting assemblies to VRML files")


def exportGLTF(assy, path: str, binary=None, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> None:
    raise NotImplementedError("exporting assemblies to glTF files")


def exportVTKJS(assy, path: str) -> None:
    raise NotImplementedError("exporting assemblies to VTK.js files")
