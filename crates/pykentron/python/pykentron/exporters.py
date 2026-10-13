"""Writing shapes to files (CadQuery's `exporters`): by type or the file's extension."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from . import _kentron as native
from .shapes import Shape


class ExportTypes:
    STL = "STL"
    STEP = "STEP"
    AMF = "AMF"
    SVG = "SVG"
    TJS = "TJS"
    DXF = "DXF"
    VRML = "VRML"
    VTP = "VTP"
    THREEMF = "3MF"
    BREP = "BREP"
    BIN = "BIN"


def _shape_of(w) -> Shape:
    """The shape of a shape, or a compound of a workplane's stack (pykentron's or
    CadQuery's) or of shapes."""
    from .shapes import Compound

    if isinstance(w, Shape):
        return w
    items = w.vals() if hasattr(w, "vals") else list(w)
    shapes = [o for o in items if isinstance(o, Shape)]
    return shapes[0] if len(shapes) == 1 else Compound.makeCompound(shapes)


def _type_of(fname, exportType) -> str:
    if exportType:
        if exportType.upper() not in vars(ExportTypes).values():
            raise ValueError("Unknown export type")
        return exportType.upper()
    suffix = Path(str(fname)).suffix.lower().lstrip(".")
    # (.binstep: STEP in binary, written as such by the native writer.)
    known = {"step": "STEP", "binstep": "STEP"}
    if suffix in known:
        return known[suffix]
    if suffix.upper() in vars(ExportTypes).values():
        return suffix.upper()
    raise ValueError("Unknown extensions, specify export type explicitly")


# Millimetres per length unit, by CadQuery's names for them.
_MM_PER = {"MM": 1.0, "CM": 10.0, "M": 1e3, "KM": 1e6, "INCH": 25.4, "IN": 25.4, "FT": 304.8, "MI": 1609344.0, "UM": 1e-3, "NM": 1e-6}


def _mm_per(unit: str) -> float:
    try:
        return _MM_PER[unit.upper()]
    except KeyError:
        raise ValueError(f"unknown length unit {unit}") from None


def export(w, fname, exportType: Optional[str] = None, tolerance: float = 0.1, angularTolerance: float = 0.1, unit: str = "MM", outputUnit=None, opt=None, **kwargs) -> None:
    """Writes a shape, or a workplane's shapes, to a file of a type (or the file's
    extension's)."""
    opt = dict(opt or {})
    shape = _shape_of(w)
    kind = _type_of(fname, exportType)
    fname = str(fname)
    if kind == "BREP":
        shape.exportBrep(fname)
    elif kind == "STEP":
        # The shape's lengths are in `unit`; the file is written in `outputUnit` (by default
        # the same), declared so.
        out = (outputUnit or unit).upper()
        factor = _mm_per(unit) / _mm_per(out)
        written = shape.scale(factor) if factor != 1.0 else shape
        native.write_step(written.wrapped, fname, Path(fname).stem, unit=out)
    elif kind == "STL":
        ascii = bool(opt.get("ascii", False) or opt.get("ASCII", False) or kwargs.get("ascii", False))
        native.write_stl(shape.wrapped, fname, tolerance, angularTolerance, ascii)
    elif kind == "SVG":
        from .occ_impl.exporters.svg import getSVG

        with open(fname, "w") as f:
            f.write(getSVG(shape, opt))
    elif kind == "TJS":
        from .occ_impl.exporters.json import JsonMesh

        vertices, triangles = shape.tessellate(tolerance, angularTolerance)
        mesher = JsonMesh()
        for v in vertices:
            mesher.addVertex(v.x, v.y, v.z)
        for t in triangles:
            mesher.addTriangleFace(*t)
        with open(fname, "w") as f:
            f.write(mesher.toJson())
    elif kind == "AMF":
        from .occ_impl.exporters.amf import AmfWriter

        writer = AmfWriter(shape.tessellate(tolerance, angularTolerance))
        with open(fname, "wb") as f:
            writer.writeAmf(f)
    elif kind == "3MF":
        from .occ_impl.exporters.threemf import ThreeMFWriter

        writer = ThreeMFWriter(shape, tolerance, angularTolerance, **opt)
        with open(fname, "wb") as f:
            writer.write3mf(f)
    elif kind == "DXF":
        from .occ_impl.exporters.dxf import exportDXF

        exportDXF(w, fname, **opt)
    else:
        raise NotImplementedError(f"exporting {kind} files")


def export_assembly(assembly, path, exportType: Optional[str] = None, tolerance: float = 0.1, angularTolerance: float = 0.1) -> None:
    kind = _type_of(path, exportType)
    if kind == "STEP":
        native.write_step_assembly(str(path), [assembly._step_tree()])
    elif kind in ("BREP", "STL"):
        export(assembly.toCompound(), path, kind, tolerance, angularTolerance)
    else:
        raise NotImplementedError(f"exporting assemblies to {kind} files")
