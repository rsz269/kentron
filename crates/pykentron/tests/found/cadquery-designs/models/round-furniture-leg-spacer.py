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
# Round furniture-leg spacer with a raised retaining rim.

import cadquery as cq
import cq_utils

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# =========================
# Parameters
# =========================
# Furniture-leg diameter inside the rim.
DIAMETER = 50.0
SPACER_HEIGHT = 5.0
RIM_SIZE = 2.0


def build_round_spacer(
    diameter: float = DIAMETER,
    spacer_height: float = SPACER_HEIGHT,
    rim_size: float = RIM_SIZE,
) -> cq.Workplane:
    """Build a round spacer sized by the furniture-leg area inside its rim."""
    if diameter <= 0.0:
        raise ValueError("Furniture-leg diameter must be positive")
    if spacer_height <= 0.0:
        raise ValueError("Spacer height must be positive")
    if rim_size <= 0.0:
        raise ValueError("Rim size must be positive")

    inner_radius = diameter / 2.0
    outer_radius = inner_radius + rim_size
    base = cq.Workplane("XY").circle(outer_radius).extrude(spacer_height)
    outer_rim = (
        cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
        .circle(outer_radius)
        .extrude(rim_size)
    )
    inner_cutout = (
        cq.Workplane("XY", origin=(0.0, 0.0, spacer_height))
        .circle(inner_radius)
        .extrude(rim_size)
    )

    return base.union(outer_rim.cut(inner_cutout))


if __name__ == "__main__":
    spacer = build_round_spacer()
    models = {
        f"Round-furniture-leg-spacer-diameter_{DIAMETER}-height_{SPACER_HEIGHT}": spacer,
    }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **models)

    # Show in OCP Viewer
    cq_utils.show_models(**models)
