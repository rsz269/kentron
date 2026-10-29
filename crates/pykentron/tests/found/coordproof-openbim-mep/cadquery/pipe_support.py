"""Parametric floor-mounted pipe support bracket."""

from __future__ import annotations

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "support_pipe_bracket_type_a"

DEFAULT_PARAMETERS = {
    "length_mm": 220,
    "width_mm": 120,
    "height_mm": 300,
    "base_thickness_mm": 16,
    "web_thickness_mm": 12,
    "pipe_diameter_mm": 60,
    "clearance_mm": 8,
    "bolt_diameter_mm": 14,
    "bolt_spacing_mm": 150,
    "material_tag": "painted_steel",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["width_mm"])
    height = float(p["height_mm"])
    base_t = float(p["base_thickness_mm"])
    web_t = float(p["web_thickness_mm"])
    pipe_d = float(p["pipe_diameter_mm"])
    clearance = float(p["clearance_mm"])
    bolt_d = float(p["bolt_diameter_mm"])
    bolt_spacing = float(p["bolt_spacing_mm"])

    base = (
        cq.Workplane("XY")
        .box(length, width, base_t)
        .faces(">Z")
        .workplane()
        .pushPoints([(-bolt_spacing / 2, 0), (bolt_spacing / 2, 0)])
        .hole(bolt_d)
    )
    post = cq.Workplane("XY").box(web_t, width * 0.55, height).translate((0, 0, height / 2))
    gusset_a = cq.Workplane("XY").box(length * 0.62, web_t, height * 0.45).translate((0, -width * 0.18, height * 0.35))
    gusset_b = cq.Workplane("XY").box(length * 0.62, web_t, height * 0.45).translate((0, width * 0.18, height * 0.35))
    saddle = cq.Workplane("XY").box(length * 0.7, width * 0.7, 34).translate((0, 0, height + 16))
    cutter = (
        cq.Workplane("YZ")
        .circle(pipe_d / 2 + clearance)
        .extrude(length)
        .translate((-length / 2, 0, height + 33))
    )
    saddle = saddle.cut(cutter)
    return base.union(post).union(gusset_a).union(gusset_b).union(saddle)


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
