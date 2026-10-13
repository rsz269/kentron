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


def _tree(assy):
    """The assembly as the native writer's nested `(name, shape, placement rows, rgb,
    children)` tuples."""
    shapes = _shapes_of(assy.obj)
    shape = None
    if shapes:
        shape = (shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)).wrapped
    rgb = tuple(assy.color.toTuple()[:3]) if assy.color is not None else None
    return (assy.name, shape, [list(r) for r in assy.loc.rows], rgb, [_tree(ch) for ch in assy.children])


def exportAssembly(assy, path: str, mode: str = "default", unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    """Writes the assembly to a STEP file, its tree as products and occurrences."""
    if mode == "fused":
        raise NotImplementedError("exporting fused assemblies")
    if unit.upper() != "MM" or (outputUnit is not None and outputUnit.upper() != "MM"):
        raise NotImplementedError("STEP units other than millimetres")
    native.write_step_assembly(str(path), [_tree(assy)])
    return True


def exportStepMeta(assy, path: str, **kwargs) -> bool:
    return exportAssembly(assy, path)


def exportCAF(assy, path: str, binary: bool = True) -> bool:
    raise NotImplementedError(f"exporting assemblies to {'XBF' if binary else 'XML'} files")


def exportVRML(assy, path: str, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> None:
    raise NotImplementedError("exporting assemblies to VRML files")


def exportGLTF(assy, path: str, binary=None, tolerance: float = 1e-3, angularTolerance: float = 0.1) -> None:
    raise NotImplementedError("exporting assemblies to glTF files")


def exportVTKJS(assy, path: str) -> None:
    raise NotImplementedError("exporting assemblies to VTK.js files")
