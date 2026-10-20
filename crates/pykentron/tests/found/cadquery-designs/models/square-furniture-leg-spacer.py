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
# Square furniture-leg spacer with a raised retaining rim that can have one or two adjacent open sides,
# if you need to fit you furniture against a wall or into a corner.

import cadquery as cq
import cq_utils

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# =========================
# Parameters
# =========================
# Furniture-leg dimensions inside the rim.
X_SIZE = 50.0
Y_SIZE = 50.0
SPACER_HEIGHT = 5.0
RIM_SIZE = 2.0
EDGE_FILLET = 0.8

# 0 keeps the complete rim, 1 opens the +Y edge for a wall, and 2 also opens
# the adjacent +X edge for a corner.
OPEN_SIDES = 0


def build_square_spacer(
    x_size: float = X_SIZE,
    y_size: float = Y_SIZE,
    spacer_height: float = SPACER_HEIGHT,
    rim_size: float = RIM_SIZE,
    edge_fillet: float = EDGE_FILLET,
    open_sides: int = OPEN_SIDES,
) -> cq.Workplane:
    """Build a rectangular spacer with zero, one, or two open rim edges."""
    if x_size <= 0.0 or y_size <= 0.0:
        raise ValueError("Furniture-leg X and Y dimensions must be positive")
    if spacer_height <= 0.0:
        raise ValueError("Spacer height must be positive")
    if rim_size <= 0.0:
        raise ValueError("Rim size must be positive")
    if edge_fillet <= 0.0 or edge_fillet >= rim_size:
        raise ValueError("Edge fillet must be positive and smaller than rim size")
    if open_sides not in (0, 1, 2):
        raise ValueError("Open sides must be 0, 1, or 2")

    outer_x_size = x_size + 2.0 * rim_size
    outer_y_size = y_size + 2.0 * rim_size
    base = cq.Workplane("XY").box(outer_x_size, outer_y_size, spacer_height, centered=(True, True, False))
    outer_rim = (
        cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
        .box(outer_x_size, outer_y_size, rim_size, centered=(True, True, False))
    )
    inner_cutout = (
        cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
        .box(x_size, y_size, rim_size, centered=(True, True, False))
    )
    rim = outer_rim.cut(inner_cutout)

    if open_sides >= 1:
        positive_y_edge = (
            cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
            .box(outer_x_size, rim_size, rim_size, centered=(True, True, False))
            .translate((0.0, (outer_y_size - rim_size) / 2.0, 0.0))
        )
        rim = rim.cut(positive_y_edge)

    if open_sides == 2:
        positive_x_edge = (
            cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
            .box(rim_size, outer_y_size, rim_size, centered=(True, True, False))
            .translate(((outer_x_size - rim_size) / 2.0, 0.0, 0.0))
        )
        rim = rim.cut(positive_x_edge)

    return base.union(rim).edges("|Z").fillet(edge_fillet)


if __name__ == "__main__":
    spacer = build_square_spacer()
    models = {
        f"Square-furniture-leg-spacer-{X_SIZE}x{Y_SIZE}-height_{SPACER_HEIGHT}-open_sides_{OPEN_SIDES}": spacer,
    }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **models)

    # Show in OCP Viewer
    cq_utils.show_models(**models)
