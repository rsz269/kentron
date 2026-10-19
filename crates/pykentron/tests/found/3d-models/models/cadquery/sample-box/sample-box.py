"""Sample box with rounded edges and a lid - demonstrates CadQuery basics."""

from pathlib import Path

import cadquery as cq

# --- Parameters ---
outer_width = 60.0  # mm - outer width of the box
outer_depth = 40.0  # mm - outer depth of the box
outer_height = 25.0  # mm - outer height of the box (excluding lid)
wall_thickness = 1.6  # mm - wall thickness (4 perimeters at 0.4mm nozzle)
corner_radius = 3.0  # mm - fillet radius for outer vertical edges
bottom_thickness = 1.2  # mm - bottom plate thickness

# --- Build ---
# Outer shell
outer = (
    cq.Workplane("XY").box(outer_width, outer_depth, outer_height).edges("|Z").fillet(corner_radius)
)

# Inner cavity (shell operation)
result = outer.faces(">Z").shell(-wall_thickness)

# --- Export ---
output_dir = Path(__file__).parent / "output"
output_dir.mkdir(exist_ok=True)
result.export(str(output_dir / "sample-box.stl"))
print(f"Exported to {output_dir / 'sample-box.stl'}")

# CQ-editor compatibility
try:
    show_object(result)  # type: ignore[name-defined]
except NameError:
    pass
