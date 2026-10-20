
# BSD 3-Clause License
#
# Copyright (c) 2025, Viktor Samokhin (wowyupiyo@gmail.com)
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
# Honeycomb infill pattern generator

import math

import cadquery as cq
from cadquery import Workplane

_AXIS_TO_PLANE = {
    "X": "YZ",
    "Y": "XZ",
    "Z": "XY",
}

_PLANE_AXES = {
    "XY": ("X", "Y"),
    "YZ": ("Y", "Z"),
    "XZ": ("X", "Z"),
}

def _hex_centers(width, height, cell_size):
    """
    Generate (x,y) centers for hex grid covering bounding box
    """
    hex_height = math.sqrt(3) * cell_size
    x_step = 1.5 * cell_size
    y_step = hex_height
    span_x = width + 2 * cell_size
    span_y = height + 2 * cell_size

    # center the lattice on origin so the pattern looks symmetric on all faces
    n_cols = max(1, int(math.ceil(span_x / x_step)) + 1)
    n_rows = max(1, int(math.ceil(span_y / y_step)) + 1)

    x0 = - (n_cols - 1) * x_step / 2.0
    y0_base = - (n_rows - 1) * y_step / 2.0

    for col in range(n_cols):
        y_offset = (hex_height / 2.0) if (col % 2 == 1) else 0.0
        y0 = y0_base + y_offset
        for row in range(n_rows):
            yield (x0 + col * x_step, y0 + row * y_step)

def _honeycomb_grid(width, height, cell_size, edge_width, thickness, plane):
    """
    Union of all hex rings covering a bounding box
    """
    centers = list(_hex_centers(width, height, cell_size))
    inner_r = max(cell_size - edge_width / 2.0, 0)
    if not centers or inner_r <= 0:
        return cq.Workplane(plane)

    pattern = (cq
               .Sketch()
               .rect(width, height)
               .push(centers)
               .regularPolygon(inner_r, 6, 90, mode="s")
    )

    return (cq
            .Workplane(plane)
            .placeSketch(pattern)
            .extrude(thickness / 2.0, both=True, combine=True)
            )

def _resolve_orientation(solid, normal_axis: str | None):
    """
    Determine which principal axis should be treated as the normal and map the
    corresponding honeycomb plane, dimensions and translation.
    """
    bb = solid.BoundingBox()  # pyright: ignore[reportAttributeAccessIssue]
    lengths = { "X": bb.xlen, "Y": bb.ylen, "Z": bb.zlen }
    centers = {
        "X": (bb.xmin + bb.xmax) / 2.0,
        "Y": (bb.ymin + bb.ymax) / 2.0,
        "Z": (bb.zmin + bb.zmax) / 2.0,
    }

    if normal_axis is None:
        normal_axis = min(lengths, key=lengths.get) # pyright: ignore[reportArgumentType, reportCallIssue]
    else:
        normal_axis = normal_axis.upper()
        if normal_axis not in _AXIS_TO_PLANE:
            raise ValueError(f"normal_axis must be one of X, Y, Z (got {normal_axis!r})")

    plane = _AXIS_TO_PLANE[normal_axis] # pyright: ignore[reportArgumentType]
    axis_u, axis_v = _PLANE_AXES[plane]

    translation = { "X": 0.0, "Y": 0.0, "Z": 0.0 }
    translation[axis_u] = centers[axis_u]
    translation[axis_v] = centers[axis_v]
    translation[normal_axis] = centers[normal_axis] # pyright: ignore[reportArgumentType]

    return {
        "plane": plane,
        "width": lengths[axis_u],
        "height": lengths[axis_v],
        "thickness": lengths[normal_axis], # pyright: ignore[reportArgumentType]
        "translation": (translation["X"], translation["Y"], translation["Z"]),
        "normal_axis": normal_axis,
    }

def apply_honeycomb(
        wall: Workplane,
        cell_size,
        edge_width,
        shell_thickness,
        normal_axis: str | None = None) -> Workplane:
    """
    Apply a honeycomb infill to any polygonal wall solid
    """
    solid = wall.val()
    orientation = _resolve_orientation(solid, normal_axis)

    # Create the outer shell following the resolved normal axis
    selector = f"+{orientation['normal_axis']} or -{orientation['normal_axis']}"
    shell_only = wall.faces(selector).shell(-shell_thickness, kind='intersection')

    # Build honeycomb grid aligned to wall orientation
    pattern = _honeycomb_grid(
        orientation["width"],
        orientation["height"],
        cell_size,
        edge_width,
        orientation["thickness"],
        orientation["plane"],
    )
    pattern = pattern.translate(orientation["translation"])

    # Intersect pattern with the fill region and wall
    infill = wall.intersect(pattern)

    return shell_only.union(infill)
