"""Parametric hollow rectangular duct segment."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "duct_main_001"

DEFAULT_PARAMETERS = {
    "length_mm": 900,
    "duct_width_mm": 420,
    "duct_height_mm": 220,
    "wall_thickness_mm": 6,
    "flange_depth_mm": 26,
    "material_tag": "galvanized_sheet_metal",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["duct_width_mm"])
    height = float(p["duct_height_mm"])
    wall = float(p["wall_thickness_mm"])
    flange = float(p["flange_depth_mm"])

    outer = cq.Workplane("XY").box(length, width, height)
    inner = cq.Workplane("XY").box(length + 4, width - 2 * wall, height - 2 * wall)
    duct = outer.cut(inner)
    flange_a = cq.Workplane("XY").box(flange, width + 36, height + 36).cut(
        cq.Workplane("XY").box(flange + 2, width, height)
    ).translate((-length / 2 - flange / 2, 0, 0))
    flange_b = cq.Workplane("XY").box(flange, width + 36, height + 36).cut(
        cq.Workplane("XY").box(flange + 2, width, height)
    ).translate((length / 2 + flange / 2, 0, 0))
    return duct.union(flange_a).union(flange_b)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
