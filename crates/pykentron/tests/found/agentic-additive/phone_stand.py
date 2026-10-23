"""
Phone Stand Generator

Generates a small, printable phone stand with:
- A flat rectangular base (stable on a table)
- An angled back support (wedge-shaped) that holds the phone

Output: phone_stand.stl
"""

import cadquery as cq
from pathlib import Path

# -------------------------------
# Dimensions (all in mm)
# -------------------------------
# Base
BASE_WIDTH = 80.0     # X dimension (matches typical phone width)
BASE_DEPTH = 60.0     # Y dimension
BASE_THICKNESS = 4.0  # Z dimension

# Back support (triangular wedge)
BACK_HEIGHT = 55.0          # How tall the support rises
BACK_DEPTH_AT_BASE = 18.0   # Thickness of the wedge at the bottom
SUPPORT_ANGLE_DEG = 70.0    # Angle of the front face from horizontal

# Phone pocket (notch where the bottom of the phone sits)
NOTCH_DEPTH = 2.5    # How deep the notch is cut
NOTCH_WIDTH = 75.0   # Width of the notch (slightly less than base)


def generate_phone_stand(output_path: Path):
    """Generate the phone stand STL."""
    # -------------------------------
    # 1. Base plate (centered on XY)
    # -------------------------------
    base = cq.Workplane("XY").box(BASE_WIDTH, BASE_DEPTH, BASE_THICKNESS)

    # -------------------------------
    # 2. Back support wedge
    # -------------------------------
    # We sketch a right-triangle profile in the YZ plane.
    # Points (y, z) relative to the sketch origin:
    #   (0, 0)               -> front bottom (meets the base)
    #   (0, BACK_HEIGHT)     -> front top (vertical face the phone rests on)
    #   (BACK_DEPTH_AT_BASE, 0) -> back bottom
    #
    # The triangle slopes backward as it rises, so the phone leans back comfortably.
    wedge = (
        cq.Workplane("YZ")
        .moveTo(0, 0)
        .lineTo(0, BACK_HEIGHT)
        .lineTo(BACK_DEPTH_AT_BASE, 0)
        .close()
        .extrude(BASE_WIDTH / 2, both=True)   # extrude +/-X to get full width
    )

    # Position the wedge so its front bottom edge sits at the back edge of the base.
    # Base extends from Y = -BASE_DEPTH/2 to +BASE_DEPTH/2.
    # We shift the wedge back by BASE_DEPTH/2 along Y, and up by BASE_THICKNESS/2 along Z
    # so it sits flush on top of the base.
    wedge = wedge.translate((0, -BASE_DEPTH / 2, BASE_THICKNESS / 2))

    # Combine base + wedge
    stand = base.union(wedge)

    # -------------------------------
    # 3. Optional notch / lip for phone bottom
    # -------------------------------
    # A small cutout right in front of the wedge gives the phone a "shelf" to sit in.
    notch = (
        cq.Workplane("XY")
        .workplane(offset=BASE_THICKNESS / 2 - 0.01)
        .center(0, -BASE_DEPTH / 2 + NOTCH_DEPTH / 2)
        .box(NOTCH_WIDTH, NOTCH_DEPTH, NOTCH_DEPTH * 2)
    )
    stand = stand.cut(notch)

    # -------------------------------
    # 4. Round outer edges for comfort / printability
    # -------------------------------
    try:
        stand = stand.edges().fillet(1.0)
    except Exception:
        pass

    # -------------------------------
    # 5. Export
    # -------------------------------
    cq.exporters.export(stand, str(output_path))
    print(f"Exported phone stand to: {output_path}")
    return output_path


if __name__ == "__main__":
    output = Path("phone_stand.stl")
    generate_phone_stand(output)
    print(f"Done! File saved as {output.absolute()}")
