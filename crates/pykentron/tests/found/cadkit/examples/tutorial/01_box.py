"""Chapter 1: a CadQuery box and loose lid. Run to write native STEP files."""
from pathlib import Path

import cadquery as cq

WIDTH = 80
DEPTH = 50
HEIGHT = 20
WALL = 2.4
LID_THICKNESS = 3


def box_body():
    # The rim is Z=0; the box extends down from it.
    outside = cq.Workplane("XY").rect(WIDTH, DEPTH).extrude(-HEIGHT)
    cavity = (cq.Workplane("XY").rect(WIDTH - 2 * WALL, DEPTH - 2 * WALL)
              .extrude(-(HEIGHT - WALL)))
    return outside.cut(cavity)


def lid_body():
    return cq.Workplane("XY").rect(WIDTH, DEPTH).extrude(LID_THICKNESS)


if __name__ == "__main__":
    output = Path("build/first-box")
    output.mkdir(parents=True, exist_ok=True)
    cq.exporters.export(box_body(), str(output / "box.step"))
    cq.exporters.export(lid_body(), str(output / "lid.step"))
    print(f"Wrote {output / 'box.step'} and {output / 'lid.step'}")
