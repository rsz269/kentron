"""
Rectangular Mounting Plate Generator

Generates a simple rectangular mounting plate, 80mm wide and 50mm tall,
with four screw holes near the corners.
"""

import cadquery as cq
from pathlib import Path

# Plate dimensions
PLATE_WIDTH = 80.0   # X dimension
PLATE_HEIGHT = 50.0  # Y dimension
PLATE_THICKNESS = 3.0  # Z dimension

# Screw hole parameters
SCREW_HOLE_DIAMETER = 4.5  # M4 screw with clearance
HOLE_INSET_X = 8.0  # Distance from left/right edges
HOLE_INSET_Y = 8.0  # Distance from top/bottom edges

# Optional: countersink parameters for flush screw heads
COUNTERSINK_DIAMETER = 8.0
COUNTERSINK_DEPTH = 1.5


def generate_mounting_plate(output_path: Path):
    """Generate the rectangular mounting plate STL."""

    # 1. Create the base plate (centered at origin)
    plate = cq.Workplane("XY").box(
        PLATE_WIDTH,
        PLATE_HEIGHT,
        PLATE_THICKNESS
    )

    # 2. Calculate hole positions (relative to center)
    # Top-right
    hx = PLATE_WIDTH / 2 - HOLE_INSET_X
    hy = PLATE_HEIGHT / 2 - HOLE_INSET_Y

    hole_positions = [
        ( hx,  hy),  # Top-right
        (-hx,  hy),  # Top-left
        ( hx, -hy),  # Bottom-right
        (-hx, -hy),  # Bottom-left
    ]

    # 3. Cut the four screw holes
    for x_pos, y_pos in hole_positions:
        hole = (
            cq.Workplane("XY")
            .workplane(offset=-PLATE_THICKNESS / 2 - 1)
            .center(x_pos, y_pos)
            .circle(SCREW_HOLE_DIAMETER / 2)
            .extrude(PLATE_THICKNESS + 2)
        )
        plate = plate.cut(hole)

    # 4. Optional: add countersinks for flush screw heads
    for x_pos, y_pos in hole_positions:
        countersink = (
            cq.Workplane("XY")
            .workplane(offset=PLATE_THICKNESS / 2 - COUNTERSINK_DEPTH)
            .center(x_pos, y_pos)
            .circle(COUNTERSINK_DIAMETER / 2)
            .extrude(COUNTERSINK_DEPTH + 0.1)
        )
        plate = plate.cut(countersink)

    # 5. Optional: round the outer edges for comfort and printability
    try:
        plate = plate.edges().fillet(1.0)
    except Exception:
        pass

    # 6. Export
    cq.exporters.export(plate, str(output_path))
    print(f"Exported mounting plate to: {output_path}")
    return output_path


if __name__ == "__main__":
    output = Path("mounting_plate.stl")
    generate_mounting_plate(output)
    print(f"Done! File saved as {output.absolute()}")
