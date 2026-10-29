"""Parametric rectangular duct trapeze hanger."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "support_duct_hanger_type_a"

DEFAULT_PARAMETERS = {
    "duct_width_mm": 420,
    "duct_height_mm": 220,
    "height_mm": 520,
    "rod_diameter_mm": 12,
    "strap_width_mm": 38,
    "strap_thickness_mm": 8,
    "material_tag": "galvanized_steel",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    duct_w = float(p["duct_width_mm"])
    duct_h = float(p["duct_height_mm"])
    height = float(p["height_mm"])
    rod_d = float(p["rod_diameter_mm"])
    strap_w = float(p["strap_width_mm"])
    strap_t = float(p["strap_thickness_mm"])
    overall_w = duct_w + 120

    top = cq.Workplane("XY").box(overall_w, strap_w, strap_t).translate((0, 0, height))
    bottom = cq.Workplane("XY").box(overall_w, strap_w, strap_t).translate((0, 0, 0))
    left_rod = cq.Workplane("XY").circle(rod_d / 2).extrude(height).translate((-overall_w / 2 + 28, 0, 0))
    right_rod = cq.Workplane("XY").circle(rod_d / 2).extrude(height).translate((overall_w / 2 - 28, 0, 0))
    duct_outline = (
        cq.Workplane("XY")
        .box(duct_w, strap_t, strap_w)
        .translate((0, -duct_h / 2, duct_h / 2))
        .union(cq.Workplane("XY").box(duct_w, strap_t, strap_w).translate((0, duct_h / 2, duct_h / 2)))
        .union(cq.Workplane("XY").box(strap_t, duct_h, strap_w).translate((-duct_w / 2, 0, duct_h / 2)))
        .union(cq.Workplane("XY").box(strap_t, duct_h, strap_w).translate((duct_w / 2, 0, duct_h / 2)))
    )
    return top.union(bottom).union(left_rod).union(right_rod).union(duct_outline)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
