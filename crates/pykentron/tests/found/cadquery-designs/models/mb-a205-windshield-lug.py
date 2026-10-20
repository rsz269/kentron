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
# v0.0.2
# Two lugs of my Mercedes-Benz wind shield (A2058680009), which is for A205 cabrio,
# were broken. There are several different types of plastic lugs in this wind shield,
# but you won't confuse them visually.
# Use a kind of reinforced with fibers filament for the print, mine was ABS-GF.
# Will see how good it holds out.

import math
from collections.abc import Sequence

import cadquery as cq
import cq_utils
from cadquery import Workplane

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# =========================
# Parameters
# =========================
BASE_W = 17.0          # X
BASE_D = 17.0          # Y
BASE_H = 20.0          # Z
BASE_FILLET_1 = 3.0
BASE_FILLET_3 = 6.0

INSERT_H = 30.0
INSERT_WALL_W = 3.25
INSERT_WALL_D = 2.5
INSERT_GROOVE_W = 6.5
INSERT_GROOVE_D = 5.5
INSERT_GROOVE_FILLET = 1.5

CUT_OFF_H = 6.0        # Z

EAR_TH = 5.0           # Thickness of ear in X
EAR_HOLE_D = 6.0       # Hole diameter
EAR_BOLT_HEAD_D = 12.5
EAR_CSK_D = 0.0
EAR_OUTER_R = 8.0      # Total length of ear from apex to the farest point of circle is 3 * R
EAR_OFFSET_X = -4.0
EAR_OFFSET_Y = 6.0

# Small rear latch/tab visible on photos (approximate)
REAR_TAB_W = 7.0
REAR_TAB_D = 3.5
REAR_TAB_H = 4.0
REAR_TAB_OFFSET_X = -2.0
REAR_TAB_OFFSET_Y = 4.0
REAR_TAB_FILLET = 0.8

RABBET_W = 2.0
RABBET_H = 4.0
RABBET_OFFSET_X1 = 6.0  # Bottom-left
RABBET_OFFSET_X2 = 2.0  # Top-right

type Point = tuple[float, float]

def create_base(
        width: float,
        depth: float,
        height: float,
        fillet_1: float,
        fillet_3: float) -> tuple[Workplane, Sequence[Point]]:

    x = width / 2.0
    y = depth / 2.0
    diff_x = width / 9.5
    p1 = (-x, y)
    p2 = (x + diff_x, y)
    p3 = (x, -y)
    p4 = (-x, -y)
    return ((
        cq.Workplane("XY")
        .polyline([p1, p2, p3, p4])
        .close()
        .extrude(height)
        .edges("|Z and >X")
        .fillet(fillet_1)
        .edges("|Z")
        .fillet(fillet_3)
    ), [p1, p2, p3, p4])

if __name__ == "__main__":
    # =========================
    # Base block
    # =========================
    part, base_p = create_base(BASE_W, BASE_D, BASE_H, BASE_FILLET_1, BASE_FILLET_3)

    insert_w = BASE_W - INSERT_WALL_W
    insert_d = BASE_D - INSERT_WALL_D
    insert, _ = create_base(
         insert_w,
         insert_d,
         INSERT_H,
         BASE_FILLET_1, BASE_FILLET_3)

    groove_left_x = -insert_w / 2.0 + INSERT_GROOVE_W / 2.0 + 2.0
    groove = (
        cq.Workplane("YZ", origin=(groove_left_x, -INSERT_H, insert_d / 2.0 - INSERT_GROOVE_D))
        .box(INSERT_H, INSERT_GROOVE_D, INSERT_GROOVE_W, centered=(False, False, True))
        .edges("<Z and |Y")
        .fillet(INSERT_GROOVE_FILLET)
        .rotate((0, 0, 0), (1, 0, 0), 90)
    )

    insert = (
        insert
        .translate((0.0, 0.0, -INSERT_H))
        .cut(groove)
        .edges("<Z or |Z")
        .fillet(0.5)
    )
    part = part.union(insert)

    # =========================
    # Ear / lug
    # =========================
    apex = (0.0, 0.0)
    left_base = (-1.5 * EAR_OUTER_R, -math.sqrt(3) * EAR_OUTER_R / 2.0)
    right_base = (-1.5 * EAR_OUTER_R, math.sqrt(3) * EAR_OUTER_R / 2.0)
    circle_center_yz = (-2.0 * EAR_OUTER_R, 0.0)

    ear_profile = (
        cq.Workplane("YZ", origin=(EAR_OFFSET_X, EAR_OFFSET_Y, BASE_H))
        .polyline([left_base, apex, right_base])
        .close()
        .pushPoints([circle_center_yz])
        .circle(EAR_OUTER_R)
        .extrude(EAR_TH)
        .faces(">X")
        .workplane()
        .pushPoints([(-2.0 * EAR_OUTER_R, 0.0)])
        #.edges("|X and (>Z or <Z)")
        #.fillet(EAR_OUTER_RADIUS - 1.0)
    )

    # =========================
    # Holes
    # =========================
    ear_profile = ear_profile.cskHole(EAR_HOLE_D, EAR_CSK_D, 82.0) if EAR_CSK_D > 0.0 else ear_profile.hole(EAR_HOLE_D)

    part = (
        part
        .faces("<X")
        .workplane()
        .pushPoints([(-EAR_OFFSET_Y + 2.0 * EAR_OUTER_R, BASE_H)])
        .hole(EAR_BOLT_HEAD_D)
        .pushPoints([(-EAR_OFFSET_Y + 2.0 * EAR_OUTER_R, BASE_H)])
        .hole(2.5 * EAR_OUTER_R, BASE_W / 2.0 + EAR_OFFSET_X)
        .union(ear_profile)
    )

    # =========================
    # L-shaped cut-off
    # =========================
    cut_off = (
        cq.Workplane("XY", origin=(0.0, 0.0, BASE_H - CUT_OFF_H))
        .polyline([base_p[0], (EAR_OFFSET_X, base_p[0][1]), (EAR_OFFSET_X, base_p[3][1]), base_p[3]])
        .close()
        .polyline([(EAR_OFFSET_X, base_p[0][1]), base_p[1], (base_p[1][0], REAR_TAB_OFFSET_Y), (EAR_OFFSET_X, REAR_TAB_OFFSET_Y)])
        .close()
        .extrude(CUT_OFF_H + 3)
    )

    part = part.cut(cut_off)

    # =========================
    # Small R-shaped tab / stop
    # =========================
    tab_p1 = (0.0, REAR_TAB_H)
    tab_p2 = (REAR_TAB_D, REAR_TAB_H)
    tab_p3 = (REAR_TAB_D, REAR_TAB_H / 2.0)
    tab_p4 = (REAR_TAB_D * 0.6, REAR_TAB_H / 2.0)
    tab_p5 = (REAR_TAB_D * 0.6, 0.0)
    tab_p6 = (0.0, 0.0)
    rear_tab = (
        cq.Workplane("YZ", origin=(REAR_TAB_OFFSET_X, REAR_TAB_OFFSET_Y, BASE_H - CUT_OFF_H))
        .polyline([tab_p1, tab_p2, tab_p3, tab_p4, tab_p5, tab_p6])
        .close()
        .extrude(REAR_TAB_W)
        .faces(">Y")
        .edges("|X ")
        .fillet(REAR_TAB_FILLET)
    )

    part = part.union(rear_tab)

    # =========================
    # Rabbet
    # =========================
    rabbet_left_x = -BASE_W / 2.0 + RABBET_OFFSET_X1 + RABBET_W / 2.0
    rabbet_right_x = BASE_W / 2.0 - RABBET_OFFSET_X2 - RABBET_W / 2.0
    rabbet_diff_x = rabbet_right_x - rabbet_left_x
    rabbet_center_x = rabbet_left_x + rabbet_diff_x / 2.0
    rabbet_p1 = (rabbet_left_x, 0.0)
    rabbet_p2 = (rabbet_left_x, BASE_H * 0.15)
    rabbet_p3 = (rabbet_center_x, BASE_H * 0.5)
    rabbet_p4 = (rabbet_right_x , BASE_H * 0.85)
    rabbet_p5 = (rabbet_right_x, BASE_H)
    rabbet_offset_z = -BASE_D / 2.0 + RABBET_W / 2.0
    rabbet_path = (
        cq.Workplane("XZ", origin=(0.0, rabbet_offset_z))
        .spline([rabbet_p1, rabbet_p2, rabbet_p3, rabbet_p4, rabbet_p5])
        #.polyline([rabbet_p1, rabbet_p2, rabbet_p3, rabbet_p4, rabbet_p5])
    )

    rabbet = (
        cq.Workplane("XY", origin=(rabbet_p1[0], rabbet_offset_z))
        .rect(RABBET_W, RABBET_H)
        .sweep(rabbet_path, isFrenet=False)
    )

    part = part.cut(rabbet)

    # =========================
    # Final cleanup
    # =========================
    #part = part.faces(">Z").edges().fillet(0.2)

    lug2 = part.mirror()

    model = { "MB a205 windshield lug 1": part, "MB a205 windshield lug 2": lug2 }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **model)

    # Show in OCP Viewer
    cq_utils.show_models(**model)
