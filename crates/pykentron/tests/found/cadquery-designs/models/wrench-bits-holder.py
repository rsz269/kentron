# BSD 3-Clause License
#
# Copyright (c) 2026, Viktor Samokhin (wowyupiyo@gmail.com)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
################################################################################
#
# v0.0.1
# Sliding drawer wrench bits holder. The drawer itself is a bit tray with
# per-bit cylindrical pockets and a rounded front cut-out for finger grip.

import math
import sys

import cadquery as cq
import cq_utils

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# =========================
# Parameters
# =========================
BIT_SPECS = [
    (22.0, 36.0),
    (18.0, 36.0),
    (22.0, 34.5),
    (16.0, 28.0),
    (14.0, 26.0),
    (17.0, 26.5),
]

# Set to 0 for automatic near-square packing.
LAYOUT_COLUMNS = 3

DRAWER_WALL_THICKNESS = 3.0
DRAWER_MIN_CELL_DEPTH = 8.0
CELL_DIAMETER_CLEARANCE = 0.6
CELL_GAP = 2.0
DRAWER_FRAME_MARGIN = 2.5
DRAWER_GRIP_DIAMETER = 22.0
DRAWER_GRIP_VERTICAL_CENTER_P = 0.62

HOUSING_WALL_THICKNESS = DRAWER_WALL_THICKNESS
HOUSING_CLEARANCE = 0.4
SLIDER_DEPTH = 3.0 * HOUSING_CLEARANCE
SLIDER_LATCH_RADIUS = 2.0

FILLET = 0.4

def resolve_columns(bit_count: int, columns: int) -> int:
    """Pick a sensible packing width if the user requests auto layout."""
    return columns if columns and columns > 0 else max(1, math.ceil(math.sqrt(bit_count)))

def fillet_external_edges(part: cq.Workplane, radius: float = FILLET) -> cq.Workplane:
    """Fillet the straight outer edges of a box-like part."""
    return part.edges("|X or |Y or |Z").fillet(radius)

def arrange_bit_specs(
    bit_specs: list[tuple[float, float]],
    columns: int,
    cell_gap: float = CELL_GAP,
    drawer_frame_margin: float = DRAWER_FRAME_MARGIN,
    drawer_min_cell_depth: float = DRAWER_MIN_CELL_DEPTH,
) -> tuple[float, float, float, float, list[tuple[float, float, float, float]]]:
    """Return layout metadata and per-bit placement info"""

    resolved_columns = resolve_columns(len(bit_specs), columns)
    rows = math.ceil(len(bit_specs) / resolved_columns)

    max_diameter = 0.0
    max_height = 0.0
    min_height = sys.float_info.max
    for diameter, height in bit_specs:
        max_diameter = max(diameter, max_diameter)
        max_height = max(height, max_height)
        min_height = min(height, min_height)

    assert min_height - 5.0 >= drawer_min_cell_depth, (
        "The shortes bit must be at least 5 mm higher than the minimum cell depth to ensure a good fit and easy removal"
    )

    pitch = max_diameter + cell_gap
    footprint_x = (resolved_columns - 1) * pitch + max_diameter
    footprint_y = (rows - 1) * pitch + max_diameter

    outer_margin = 2.0 * drawer_frame_margin
    drawer_outer_x = footprint_x + outer_margin
    drawer_outer_y = footprint_y + outer_margin

    x_origin = -((resolved_columns - 1) * pitch) / 2.0
    y_origin = -((rows - 1) * pitch) / 2.0

    placements = []
    for index, (diameter, height) in enumerate(bit_specs):
        row = index // resolved_columns
        col = index % resolved_columns
        placements.append((
            x_origin + col * pitch,
            y_origin + row * pitch,
            diameter,
            height
        ))

    return (
        min_height,
        max_height,
        drawer_outer_x,
        drawer_outer_y,
        placements,
    )

def build_drawer(
    bit_specs: list[tuple[float, float]],
    columns: int,
    drawer_wall_thickness: float = DRAWER_WALL_THICKNESS,
    drawer_min_cell_depth: float = DRAWER_MIN_CELL_DEPTH,
    pocket_diameter_clearance: float = CELL_DIAMETER_CLEARANCE,
    drawer_grip_diameter: float = DRAWER_GRIP_DIAMETER,
    housing_wall_thickness: float = HOUSING_WALL_THICKNESS,
) -> tuple[cq.Workplane, float, float, float, float]:
    """Build the sliding drawer that also acts as the bit tray."""
    (
        min_height,
        max_height,
        drawer_outer_x,
        drawer_outer_y,
        placements,
    ) = arrange_bit_specs(bit_specs, columns)

    drawer_outer_z = max_height + drawer_wall_thickness
    drawer_inner_z = max_height - min_height + drawer_min_cell_depth + drawer_wall_thickness
    drawer = (
        cq
        .Workplane("XY", origin=(0.0, 0.0, housing_wall_thickness))
        .box(drawer_outer_x, drawer_outer_y, drawer_inner_z, centered=(True, True, False))
        )

    delta_z = drawer_outer_z - drawer_inner_z
    assert delta_z - drawer_grip_diameter / 2.0 >= 5.0, (
        "The grip hole must be at least 5 mm away from the inner face of the drawer to ensure structural integrity"
    )
    for x, y, diameter, height in placements:
         drawer = (
             drawer
             .faces(">Z")
             .workplane(origin=(x, y))
             .hole(diameter + pocket_diameter_clearance, height - delta_z)
         )

    front = (
        cq
        .Workplane("YZ", origin=(-drawer_outer_x / 2.0 - drawer_wall_thickness, 0.0, housing_wall_thickness))
        .box(drawer_outer_y, drawer_outer_z, drawer_wall_thickness, centered=(True, False, False))
        )

    grip_center_z = drawer_outer_z + housing_wall_thickness
    front = (
        front
        .faces("<X")
        .workplane(origin=(-drawer_outer_x / 2.0 + drawer_wall_thickness / 2.0, 0.0, grip_center_z))
        .hole(drawer_grip_diameter, drawer_wall_thickness)
    )

    return drawer.union(front), drawer_outer_x, drawer_outer_y, drawer_outer_z, drawer_inner_z

def add_sliders(
        housing: cq.Workplane,
        drawer: cq.Workplane,
        drawer_outer_x: float,
        drawer_outer_y: float,
        drawer_inner_z: float,
        housing_wall_thickness: float = HOUSING_WALL_THICKNESS,
        drawer_wall_thickness: float = DRAWER_WALL_THICKNESS,
        slider_depth: float = SLIDER_DEPTH,
        housing_clearance: float = HOUSING_CLEARANCE,
        sldier_latch_radius: float = SLIDER_LATCH_RADIUS) -> tuple[cq.Workplane, cq.Workplane]:
    workplane_x = -drawer_outer_x / 2.0 - drawer_wall_thickness
    workplane_y = drawer_outer_y / 2.0
    latch_origin = (-drawer_outer_x / 2.0, 0.0, drawer_inner_z / 2.0 + housing_wall_thickness)
    latch_offset = -sldier_latch_radius * 0.75
    left_slider = (
        cq
        .Workplane("XY", origin=(workplane_x, workplane_y, housing_wall_thickness))
        .box(drawer_outer_x + drawer_wall_thickness, slider_depth, drawer_inner_z, centered=(False, False, False))
        .faces(">Y")
        .workplane(offset=latch_offset, origin=latch_origin)
        .sphere(sldier_latch_radius)
        )
    drawer = drawer.union(left_slider).union(left_slider.mirror("XZ"))

    left_cutout = (
        cq
        .Workplane("XY", origin=(workplane_x, workplane_y, housing_wall_thickness))
        .box(drawer_outer_x + drawer_wall_thickness, slider_depth + housing_clearance, drawer_inner_z + housing_clearance, centered=(False, False, False))
                .faces(">Y")
        .workplane(offset=latch_offset, origin=latch_origin)
        .sphere(sldier_latch_radius)
        )
    housing = housing.cut(left_cutout).cut(left_cutout.mirror("XZ"))

    return housing, drawer

def build_housing(
    drawer_outer_x: float,
    drawer_outer_y: float,
    drawer_outer_z: float,
    housing_clearance: float = HOUSING_CLEARANCE,
    housing_wall_thickness: float = HOUSING_WALL_THICKNESS,
) -> cq.Workplane:
    """Build a box-like sleeve that the drawer slides in."""
    housing_inner_x = drawer_outer_x
    housing_inner_y = drawer_outer_y + 2.0 * housing_clearance
    housing_inner_z = drawer_outer_z + housing_clearance

    walls_margin = 2.0 * housing_wall_thickness
    housing_outer_x = housing_inner_x + walls_margin
    housing_outer_y = housing_inner_y + walls_margin
    housing_outer_z = housing_inner_z + walls_margin

    outer = cq.Workplane("XY").box(housing_outer_x, housing_outer_y, housing_outer_z, centered=(True, True, False))
    outer = fillet_external_edges(outer)

    inner = (
        cq.Workplane("XY")
        .box(housing_inner_x + housing_wall_thickness, housing_inner_y, housing_inner_z, centered=(True, True, False))
        .translate((-housing_wall_thickness / 2.0, 0.0, housing_wall_thickness))
    )

    return outer.cut(inner)

if __name__ == "__main__":
    drawer, drawer_outer_x, drawer_outer_y, drawer_outer_z, drawer_inner_z = build_drawer(
        BIT_SPECS,
        LAYOUT_COLUMNS,
    )
    housing = build_housing(drawer_outer_x, drawer_outer_y, drawer_outer_z)
    housing, drawer = add_sliders(housing, drawer, drawer_outer_x, drawer_outer_y, drawer_inner_z)

    drawer = fillet_external_edges(drawer)

    model_name = (
        f"wrench-bits-holder-drawer-bits_{len(BIT_SPECS)}"
        f"-drawer_l_{drawer_outer_x}-drawer_w_{drawer_outer_y}-drawer_h_{drawer_inner_z}"
    )

    models = {
        f"{model_name}-housing": housing,
        f"{model_name}-drawer": drawer
    }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **models)

    # Show in OCP Viewer
    cq_utils.show_models(**models)
