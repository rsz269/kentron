"""Parametric wall penetration sleeve."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "sleeve_wall_penetration_type_a"

DEFAULT_PARAMETERS = {
    "pipe_diameter_mm": 110,
    "wall_thickness_mm": 200,
    "clearance_mm": 30,
    "sleeve_wall_mm": 6,
    "flange_diameter_mm": 190,
    "flange_thickness_mm": 8,
    "material_tag": "steel_sleeve",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    pipe_d = float(p["pipe_diameter_mm"])
    wall_t = float(p["wall_thickness_mm"])
    clearance = float(p["clearance_mm"])
    sleeve_wall = float(p["sleeve_wall_mm"])
    flange_d = float(p["flange_diameter_mm"])
    flange_t = float(p["flange_thickness_mm"])
    inner_d = pipe_d + clearance
    outer_d = inner_d + sleeve_wall * 2

    sleeve = cq.Workplane("XY").circle(outer_d / 2).circle(inner_d / 2).extrude(wall_t)
    flange_a = cq.Workplane("XY").circle(flange_d / 2).circle(inner_d / 2).extrude(flange_t)
    flange_b = cq.Workplane("XY").circle(flange_d / 2).circle(inner_d / 2).extrude(flange_t).translate((0, 0, wall_t - flange_t))
    return sleeve.union(flange_a).union(flange_b)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
