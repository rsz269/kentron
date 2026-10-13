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


def _tree(assy, factor: float = 1.0):
    """The assembly as the native writer's nested `(name, shape, placement rows, rgb,
    children)` tuples, lengths scaled by `factor` (shapes and placements' translations)."""
    shapes = _shapes_of(assy.obj)
    shape = None
    if shapes:
        whole = shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)
        shape = (whole.scale(factor) if factor != 1.0 else whole).wrapped
    rgb = tuple(assy.color.toTuple()[:3]) if assy.color is not None else None
    rows = [list(r[:3]) + [r[3] * factor] for r in assy.loc.rows]
    return (assy.name, shape, rows, rgb, [_tree(ch, factor) for ch in assy.children])


def exportAssembly(assy, path: str, mode: str = "default", unit: str = "MM", outputUnit=None, **kwargs) -> bool:
    """Writes the assembly to a STEP file, its tree as products and occurrences."""
    if mode == "fused":
        raise NotImplementedError("exporting fused assemblies")
    from ...exporters import _mm_per

    out = (outputUnit or unit).upper()
    native.write_step_assembly(str(path), [_tree(assy, _mm_per(unit) / _mm_per(out))], unit=out)
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
