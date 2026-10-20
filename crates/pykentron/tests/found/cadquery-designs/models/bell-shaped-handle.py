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
# Bell-shaped hanging handle? Does it tell you something? :-)
# It's a thing you can put a lace through its top hole, tie a knot inside the handle
# and you have it.

import cadquery as cq
import cq_utils

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# =========================
# Parameters
# =========================
HEIGHT = 50.0

BOTTOM_D = 40.0
TOP_D = 15.0
TOP_HOLE_D = 5.0

SHOULDER_HEIGHT_P = 0.82
SHOULDER_DIAMETER_P = 0.62

WALL_TH = 3.0
FILLET = 0.8

if __name__ == "__main__":
    bottom_r = BOTTOM_D / 2.0
    shoulder_r = BOTTOM_D * SHOULDER_DIAMETER_P / 2.0
    top_r = TOP_D / 2.0

    # Revolve a half-section around the Z axis to get a smooth bell silhouette.
    outer = (
        cq.Workplane("XZ")
        .moveTo(0.0, 0.0)
        .lineTo(bottom_r, 0.0)
        .spline(
            [
                (bottom_r, 0.0),
                (shoulder_r, HEIGHT * SHOULDER_HEIGHT_P),
                (top_r, HEIGHT),
            ],
            #tangents=((1.0, 0.1), (-0.4, 1.0)),
        )
        .lineTo(0.0, HEIGHT)
        .close()
        .revolve(360.0, (0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    )

    handle = (
        outer
        .faces("<Z")
        .shell(-WALL_TH)
        .faces(">Z")
        .workplane()
        .hole(TOP_HOLE_D)
        .edges("<Z or >Z")
        .fillet(FILLET)
    )

    model = { f'Bell-shaped hanging handle-height_{HEIGHT}-bottom_d_{BOTTOM_D}-top_d_{TOP_D}': handle }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **model)

    # Show in OCP Viewer
    cq_utils.show_models(**model)
