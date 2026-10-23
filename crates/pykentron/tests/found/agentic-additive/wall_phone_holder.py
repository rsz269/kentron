"""
Wall-Mounted Phone Holder Generator

Generates a wall-mounted phone holder for a phone approximately
75mm wide and 10mm thick, with some clearance.
"""

import cadquery as cq
from pathlib import Path

# Phone dimensions
PHONE_WIDTH = 75.0
PHONE_THICKNESS = 10.0

# Clearances
WIDTH_CLEARANCE = 2.0   # Total: 77mm slot
THICKNESS_CLEARANCE = 2.0  # Total: 12mm slot

# Holder dimensions
SLOT_WIDTH = PHONE_WIDTH + WIDTH_CLEARANCE        # 77 mm
SLOT_DEPTH = PHONE_THICKNESS + THICKNESS_CLEARANCE # 12 mm
SLOT_HEIGHT = 50.0                                 # How deep the pocket is
WALL_THICKNESS = 3.0                               # Wall thickness
BACKPLATE_WIDTH = SLOT_WIDTH + WALL_THICKNESS * 2  # 83 mm
BACKPLATE_HEIGHT = SLOT_HEIGHT + WALL_THICKNESS * 2 + 30.0  # Extra height above pocket
BACKPLATE_THICKNESS = 3.0

# Mounting holes
MOUNT_HOLE_DIAMETER = 4.5
MOUNT_HOLE_SPACING_V = 50.0  # Vertical spacing between holes
MOUNT_HOLE_INSET = 10.0      # Distance from top/bottom edges

# Bottom lip to prevent phone from sliding out
LIP_HEIGHT = 5.0
LIP_THICKNESS = 2.0

# Fillet radius
FILLET_R = 2.0


def generate_wall_phone_holder(output_path: Path):
    """Generate the wall-mounted phone holder STL."""

    # 1. Main backplate
    backplate = cq.Workplane("XY").box(
        BACKPLATE_WIDTH,
        BACKPLATE_HEIGHT,
        BACKPLATE_THICKNESS
    )

    # 2. Create the U-shaped pocket (open at top)
    # The pocket sits on the front face of the backplate
    pocket_outer = (
        cq.Workplane("XY")
        .workplane(offset=BACKPLATE_THICKNESS)
        .box(SLOT_WIDTH + WALL_THICKNESS * 2, SLOT_HEIGHT + WALL_THICKNESS, SLOT_DEPTH + WALL_THICKNESS)
    )

    # Cut out the inner slot (open at top)
    # The slot goes from the bottom of the pocket up to the top
    slot = (
        cq.Workplane("XY")
        .workplane(offset=BACKPLATE_THICKNESS + WALL_THICKNESS)
        .box(SLOT_WIDTH, SLOT_HEIGHT + 1, SLOT_DEPTH + 1)
    )

    pocket = pocket_outer.cut(slot)

    # 3. Combine backplate and pocket
    holder = backplate.union(pocket)

    # 4. Add a small bottom lip inside the pocket to keep phone from sliding out
    lip = (
        cq.Workplane("XY")
        .workplane(offset=BACKPLATE_THICKNESS + WALL_THICKNESS)
        .box(SLOT_WIDTH, LIP_HEIGHT, LIP_THICKNESS)
        .translate((0, -SLOT_HEIGHT / 2 + LIP_HEIGHT / 2, SLOT_DEPTH / 2 - LIP_THICKNESS / 2))
    )
    holder = holder.union(lip)

    # 5. Mounting holes (countersunk or through-holes)
    # Position: centered horizontally, spaced vertically near top and bottom
    hole_z_top = BACKPLATE_HEIGHT / 2 - MOUNT_HOLE_INSET
    hole_z_bottom = -BACKPLATE_HEIGHT / 2 + MOUNT_HOLE_INSET

    for hz in [hole_z_top, hole_z_bottom]:
        hole = (
            cq.Workplane("YZ")
            .workplane(offset=0)
            .center(hz, 0)
            .circle(MOUNT_HOLE_DIAMETER / 2)
            .extrude(BACKPLATE_THICKNESS + 2, both=True)
        )
        holder = holder.cut(hole)

    # 6. Round edges for comfort and printability
    try:
        holder = holder.edges("|X").fillet(FILLET_R)
    except Exception:
        pass
    try:
        holder = holder.edges("|Y").fillet(FILLET_R)
    except Exception:
        pass

    # 7. Export
    cq.exporters.export(holder, str(output_path))
    print(f"Exported wall-mounted phone holder to: {output_path}")
    return output_path


if __name__ == "__main__":
    output = Path("wall_phone_holder.stl")
    generate_wall_phone_holder(output)
    print(f"Done! File saved as {output.absolute()}")
