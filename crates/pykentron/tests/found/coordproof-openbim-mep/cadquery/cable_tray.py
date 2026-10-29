"""Parametric perforated cable tray segment."""

from __future__ import annotations

import math

from asset_io import export_shape, merged

import cadquery as cq

ASSET_ID = "cable_tray_overhead_001"

DEFAULT_PARAMETERS = {
    "length_mm": 900,
    "width_mm": 220,
    "height_mm": 80,
    "thickness_mm": 4,
    "hole_diameter_mm": 20,
    "hole_pitch_mm": 110,
    "material_tag": "perforated_aluminum",
}


def build(parameters: dict[str, object] | None = None) -> cq.Workplane:
    p = merged(parameters, DEFAULT_PARAMETERS)
    length = float(p["length_mm"])
    width = float(p["width_mm"])
    height = float(p["height_mm"])
    thickness = float(p["thickness_mm"])
    hole_d = float(p["hole_diameter_mm"])
    pitch = float(p["hole_pitch_mm"])

    x_positions = [x for x in _frange(-length / 2 + pitch, length / 2 - pitch, pitch)]
    base = (
        cq.Workplane("XY")
        .box(length, width, thickness)
        .faces(">Z")
        .workplane()
        .pushPoints([(x, -width * 0.25) for x in x_positions] + [(x, width * 0.25) for x in x_positions])
        .hole(hole_d)
    )
    left = cq.Workplane("XY").box(length, thickness, height).translate((0, -width / 2 + thickness / 2, height / 2))
    right = cq.Workplane("XY").box(length, thickness, height).translate((0, width / 2 - thickness / 2, height / 2))
    return base.union(left).union(right)


def _frange(start: float, stop: float, step: float) -> list[float]:
    if not math.isfinite(step) or step <= 0:
        raise ValueError("step must be a finite positive number")
    if not math.isfinite(start) or not math.isfinite(stop):
        raise ValueError("start and stop must be finite")
    values: list[float] = []
    value = start
    while value <= stop:
        values.append(value)
        value += step
    return values


if __name__ == "__main__":
    export_shape(build(), ASSET_ID)
