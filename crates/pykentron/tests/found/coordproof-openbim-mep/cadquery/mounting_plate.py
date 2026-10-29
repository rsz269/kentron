"""Parametric mounting plate with slotted holes."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "plate_mounting_type_a"

DEFAULT_PARAMETERS = {
    "length_mm": 260,
    "width_mm": 160,
    "thickness_mm": 12,
    "bolt_diameter_mm": 14,
    "bolt_spacing_mm": 190,
    "slot_length_mm": 34,
    "material_tag": "steel",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["width_mm"])
    thickness = float(p["thickness_mm"])
    bolt_d = float(p["bolt_diameter_mm"])
    bolt_spacing = float(p["bolt_spacing_mm"])
    slot_length = float(p["slot_length_mm"])

    numeric = {
        "length_mm": length,
        "width_mm": width,
        "thickness_mm": thickness,
        "bolt_diameter_mm": bolt_d,
        "bolt_spacing_mm": bolt_spacing,
        "slot_length_mm": slot_length,
    }
    invalid = [name for name, value in numeric.items() if value <= 0]
    if invalid:
        raise ValueError(f"Mounting plate dimensions must be positive: {', '.join(invalid)}")
    if slot_length < bolt_d:
        raise ValueError("slot_length_mm must be at least bolt_diameter_mm")
    if bolt_spacing + slot_length >= length:
        raise ValueError("bolt spacing and slot length must fit inside the plate length")

    plate = cq.Workplane("XY").box(length, width, thickness).edges("|Z").chamfer(5)
    cutter = (
        cq.Workplane("XY")
        .pushPoints([(-bolt_spacing / 2, 0), (bolt_spacing / 2, 0)])
        .slot2D(slot_length, bolt_d, 0)
        .extrude(thickness + 4)
        .translate((0, 0, -thickness / 2 - 2))
    )
    return plate.cut(cutter)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
