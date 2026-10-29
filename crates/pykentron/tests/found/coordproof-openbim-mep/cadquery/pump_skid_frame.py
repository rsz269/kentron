"""Parametric pump skid steel frame."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "equipment_pump_skid_001"

DEFAULT_PARAMETERS = {
    "length_mm": 1100,
    "width_mm": 520,
    "height_mm": 180,
    "rail_width_mm": 70,
    "crossmember_count": 3,
    "material_tag": "painted_steel",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["width_mm"])
    height = float(p["height_mm"])
    rail = float(p["rail_width_mm"])
    count = int(p["crossmember_count"])

    left_rail = cq.Workplane("XY").box(length, rail, height).translate((0, -width / 2 + rail / 2, height / 2))
    right_rail = cq.Workplane("XY").box(length, rail, height).translate((0, width / 2 - rail / 2, height / 2))
    frame = left_rail.union(right_rail)
    if count > 1:
        pitch = length / (count - 1)
        positions = [-length / 2 + i * pitch for i in range(count)]
    else:
        positions = [0]
    for x in positions:
        frame = frame.union(cq.Workplane("XY").box(rail, width, height * 0.55).translate((x, 0, height * 0.275)))
    return frame.edges("|Z").chamfer(4)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
