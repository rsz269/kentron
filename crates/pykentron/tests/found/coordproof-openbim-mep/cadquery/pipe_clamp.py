"""Parametric pipe clamp with mounting ears."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "pipe_clamp_type_a"

DEFAULT_PARAMETERS = {
    "pipe_diameter_mm": 60,
    "width_mm": 45,
    "thickness_mm": 10,
    "ear_length_mm": 34,
    "bolt_diameter_mm": 10,
    "material_tag": "galvanized_steel",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    pipe_d = float(p["pipe_diameter_mm"])
    width = float(p["width_mm"])
    thickness = float(p["thickness_mm"])
    ear_length = float(p["ear_length_mm"])
    bolt_d = float(p["bolt_diameter_mm"])
    inner_r = pipe_d / 2
    outer_r = inner_r + thickness

    ring = cq.Workplane("XY").circle(outer_r).circle(inner_r).extrude(width)
    left_ear = cq.Workplane("XY").box(ear_length, thickness * 2.2, width).translate((-outer_r - ear_length / 2 + 2, 0, width / 2))
    right_ear = cq.Workplane("XY").box(ear_length, thickness * 2.2, width).translate((outer_r + ear_length / 2 - 2, 0, width / 2))
    clamp = ring.union(left_ear).union(right_ear)
    bolt_cutters = (
        cq.Workplane("XY")
        .pushPoints([(-outer_r - ear_length / 2 + 2, 0), (outer_r + ear_length / 2 - 2, 0)])
        .circle(bolt_d / 2)
        .extrude(width + 4)
        .translate((0, 0, -2))
    )
    return clamp.cut(bolt_cutters)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
