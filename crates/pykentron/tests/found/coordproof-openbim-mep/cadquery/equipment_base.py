"""Parametric concrete equipment housekeeping pad."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "equipment_base_type_a"

DEFAULT_PARAMETERS = {
    "length_mm": 1200,
    "width_mm": 800,
    "height_mm": 160,
    "bolt_diameter_mm": 18,
    "bolt_offset_mm": 120,
    "material_tag": "concrete",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["width_mm"])
    height = float(p["height_mm"])
    bolt_d = float(p["bolt_diameter_mm"])
    offset = float(p["bolt_offset_mm"])
    points = [
        (-length / 2 + offset, -width / 2 + offset),
        (length / 2 - offset, -width / 2 + offset),
        (-length / 2 + offset, width / 2 - offset),
        (length / 2 - offset, width / 2 - offset),
    ]
    return (
        cq.Workplane("XY")
        .box(length, width, height)
        .edges("|Z")
        .chamfer(12)
        .faces(">Z")
        .workplane()
        .pushPoints(points)
        .hole(bolt_d, depth=height * 0.65)
    )


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
