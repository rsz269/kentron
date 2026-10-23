from __future__ import annotations

from pathlib import Path
import sys

import cadquery as cq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from cupboard_benchmark.exports import render_preview_png, write_manifest, write_obj
from cupboard_benchmark.spec import cupboard_parts


OUT = ROOT / "outputs" / "cadquery"


def make_box(part):
    cx, cy, cz = part.center
    return cq.Workplane("XY").box(part.dx, part.dy, part.dz).translate((cx, cy, cz))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    parts = cupboard_parts()
    solids = []
    assembly = cq.Assembly(name="cadquery_cupboard")
    for part in parts:
        shape = make_box(part)
        solids.append(shape.val())
        r, g, b = [channel / 255 for channel in part.color]
        assembly.add(shape, name=part.name, color=cq.Color(r, g, b))

    compound = cq.Compound.makeCompound(solids)
    step_path = OUT / "cupboard.step"
    stl_path = OUT / "cupboard.stl"
    obj_path = OUT / "cupboard.obj"
    preview_path = OUT / "cupboard_preview.png"

    try:
        assembly.save(str(step_path))
    except Exception:
        cq.exporters.export(compound, str(step_path))
    cq.exporters.export(compound, str(stl_path), tolerance=0.2, angularTolerance=0.2)
    write_obj(parts, obj_path)
    render_preview_png(parts, preview_path, "CadQuery cupboard benchmark")

    write_manifest(
        "CadQuery",
        OUT,
        {
            "source": "scripts/generate_cadquery.py",
            "step": str(step_path.relative_to(ROOT)),
            "stl": str(stl_path.relative_to(ROOT)),
            "obj": str(obj_path.relative_to(ROOT)),
            "preview": str(preview_path.relative_to(ROOT)),
        },
        ["CadQuery exports native STEP and STL through OpenCascade/OCP."],
    )


if __name__ == "__main__":
    main()
