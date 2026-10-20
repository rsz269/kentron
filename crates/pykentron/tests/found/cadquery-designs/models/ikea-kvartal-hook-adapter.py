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
# Hook adapter for the IKEA Kvartal gliders.
# Print it slowing down and intensifying cooling on overhangs.

import cadquery as cq
import cq_utils

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

if __name__ == "__main__":
    base_height = 1.25
    base_diameter = 7.2
    base_fillet = 0.1

    mid_height = 2.3
    mid_diameter = 4.9

    ring_plate_height = 4.5
    ring_base_x = 13
    ring_base_y = 5

    torus_id = 6
    torus_od = 10
    torus_embed = 3  # how much the ring sinks into the base

    tube_radius = (torus_od - torus_id) / 4  # half the difference, divided again for radius
    major_radius = torus_od / 2 - tube_radius

    ring_plate_top_height = ring_plate_height / 4.0

    hook = (
        cq.Workplane("XY")
        .circle(base_diameter / 2)
        .extrude(base_height)
        .edges("not |Z")
        .fillet(base_fillet)
        .faces(">Z")
        .workplane()
        .circle(mid_diameter / 2)
        .extrude(mid_height)
        .faces(">Z")
        .workplane()
        .circle(mid_diameter / 2)
        .workplane(offset=ring_plate_height - ring_plate_top_height)
        .rect(ring_base_x, ring_base_y)
        .loft()
        .faces(">Z")
        .workplane()
        .rect(ring_base_x, ring_base_y)
        .extrude(ring_plate_top_height)
        .edges("|Z or >Z")
        .fillet(base_fillet)
    )

    total_height = base_height + mid_height + ring_plate_height
    torus_center_z = total_height - torus_embed + torus_od / 2

    # Vertical torus aligned to the long side (X axis)
    torus = (
        cq.Workplane(inPlane = "XY", origin = (-major_radius, 0, torus_center_z))
        .circle(tube_radius)
        .revolve(360, (major_radius, 0, 0), (major_radius, 1, 0))
    )

    hook = hook.union(torus)

    model = { "ikea_kvartal_adapter": hook}

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **model)

    # Show in OCP Viewer
    cq_utils.show_models(**model)
