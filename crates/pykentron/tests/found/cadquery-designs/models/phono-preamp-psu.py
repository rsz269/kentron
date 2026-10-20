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
# v0.0.5
# DIY phono-preamp and separate PSU cases.
# In my setup the PSU unit is only a 15-0-15V transformer, becasue rectifier circuit is built into the phono-preamp.
# For a better noise protection, the phono-preamp case should be shielded from inside using an adhesive copper foil,
# connected finally to the ground terminal.
# BOM:
# - 8x M3x8 screws + melt insert nuts for lid attachment to base
# - 8x M3 screws + nuts + PCB spacers for PCB mounting
# - 8x self-adhesive rubber feet
# - 4x cinch/RCA panel mount connectors
# - 2x 3-pole power panel mount connectors
# - 1x panel mount ground terminal
# - 1x AC power panel mount socket with fuse and switch
# - 1x power cord with plug
# - 1x MOFI-T-Phonum 2 MM/MC DIY phono-preamp kit
# - 1x suitable 15-0-15V 15VA PCB mounting transformer (e.g. Talema 70053K/72453K)
# - 1x PCB for the transformer
# - wires, cables, etc.

import cadquery as cq
import cq_utils
from cadquery import Workplane

# Toggle exports (set True when you want STL/STEP files written)
DO_STL_EXPORT = True
DO_STEP_EXPORT = False

# Configuration (tweak if needed)
PREAMP_WIDTH, PREAMP_LENGTH = 110.0, 160.0
PSU_WIDTH, PSU_LENGTH = 90.0, 135.0
HEIGHT = 50.0

BASE_WALL_THICK = 3.0              # Generic wall thickness
BASE_BOTTOM_THICK = 3.0            # Bottom plate thickness
LID_TOP_THICK = BASE_BOTTOM_THICK  # Lid top plate thickness
LID_WALL_THICK = BASE_WALL_THICK   # Side wall thickness attached to lid

PCB_SCREW_HOLE_DIAMETER = 3.2      # PCB screw should fit through
PCB_SCREW_CBR_DIAMETER = 6.5       # PCB screw counterbore diameter
PCB_SCREW_CBR_DEPTH = 2.5          # PCB screw counterbore depth
PCB_STAGE_HEIGHT = 3.1             # Stage for PCB mounts
PCB_STAGE_DIAMETER = PCB_SCREW_HOLE_DIAMETER * 1.5

LID_SCREW_HOLE_DIAMETER = 3.2      # Lid mounting screw should fit through
LID_SCREW_CBR_DIAMETER = 6.4       # Lid screw counterbore diameter
LID_SCREW_CBR_DEPTH = 1.3          # Lid screw counterbore depth

PREAMP_PCB_SPACING_X = 90.0   # Phono-preamp PCB mounting holes spacing in X direction
PREAMP_PCB_SPACING_Y = 110.0  # Phono-preamp PCB mounting holes spacing in Y direction

PSU_PCB_SPACING_X = 66.0   # PSU PCB mounting holes spacing in X direction
PSU_PCB_SPACING_Y = 71.0   # PSU PCB mounting holes spacing in Y direction

CINCH_DIAMETER = 8.2  # Diameter of cinch/RCA holes
POWER_DIAMETER = 12.5 # Diameter of power connector hole
GND_DIAMETER = 7.1    # Diameter of ground connector hole
POWER_SOCKET_WIDTH = 47.0
POWER_SOCKET_HEIGHT = 27.5

NUT_HOLE_DIAMETER = 4.0                  # Eembeded nut outer diameter
BOSS_DIAMETER = NUT_HOLE_DIAMETER * 2.0  # Attach-lid-to-base boss diameter

VENT_SPACING = 10.0    # Spacing between vents
VENT_SLOT_WIDTH = 3.0  # Width of each vent slot

FILLET_RADIUS = 0.5      # General fillet radius for edges

def add_vertical_wall_vents(wall: Workplane,
                            area_length, area_height,
                            wall_thick,
                            spacing=VENT_SPACING,
                            slot_w=VENT_SLOT_WIDTH,
                            slot_h=None) -> Workplane:
    """
    Add vertical vents as rectangular slots on a vertical wall workplane
    """
    if slot_h is None:
        slot_h = max(4.0, area_height - 4.0)

    n = max(1, int(area_length // spacing))
    actual_spacing = area_length / n
    start = -area_length / 2.0 + actual_spacing / 2.0
    ys = [start + i * actual_spacing for i in range(n)]
    pts = [(y, 0) for y in ys]
    # create slots in the wall (workplane coords map to local face Y,Z)
    return (wall
            .pushPoints(pts)
            .rect(slot_w, slot_h)
            .cutBlind(-wall_thick))

def add_diagonal_wall_vents(wall: Workplane,
                            area_length, area_height,
                            wall_thick, width,
                            height=HEIGHT,
                            base_thick=BASE_BOTTOM_THICK,
                            spacing=VENT_SPACING,
                            slot_w=VENT_SLOT_WIDTH,
                            slot_h=None,
                            angle=45.0) -> Workplane:
    """
    Add diagonal vents as rectangular slots on a vertical wall workplane.
    I'm lazy to find a better solution, so for to get proper slots in >X faces (symmetrical to those in <X faces),
    simply change the sign of wall_thick, width, and angle values.
    """
    if slot_h is None:
        slot_h = area_height - 4.0

    n = max(1, int(area_length // spacing))
    actual_spacing = area_length / n
    start = -area_length / 2.0 + actual_spacing / 2.0

    x_offset = width / 2.0 - wall_thick / 2.0
    z_offset = height / 2.0 - area_height / 2.0 + base_thick
    cuts = cq.Workplane("XY")
    for i in range(n):
        y = start + i * actual_spacing
        slot = (cq
                .Workplane("XY")
                .rect(wall_thick * 2.0, slot_w)   # long thin rectangle
                .extrude(slot_h)                  # make into a thin plate
                .rotate((0, 0, 0), (0, 0, 1), angle)
                .translate((x_offset, y, z_offset)))
        cuts = cuts.union(slot)

    # subtract the slots
    return wall.cut(cuts)

def add_pcb_mounts(bottom_plate: Workplane,
                   spacing_x, spacing_y, offset_y,
                   stage_height=PCB_STAGE_HEIGHT,
                   stage_dia=PCB_STAGE_DIAMETER,
                   bottom_thick=BASE_BOTTOM_THICK,
                   thru_dia=PCB_SCREW_HOLE_DIAMETER,
                   cbr_dia=PCB_SCREW_CBR_DIAMETER,
                   cbr_depth=PCB_SCREW_CBR_DEPTH) -> Workplane:
    """
    Add PCB mounting holes in the bottom plate and extruded stages for PCB support
    """
    hx = spacing_x / 2.0
    hy = spacing_y / 2.0
    pts = [(hx, hy), (-hx, hy), (hx, -hy), (-hx, -hy)]
    stages = cq.Workplane("XY")
    for (x, y) in pts:
        stage = (cq
                 .Workplane("XY")
                 .transformed(offset=(x, y + offset_y, bottom_thick))
                 .circle(stage_dia)
                 .extrude(stage_height))
        stages = stages.union(stage)

    return (bottom_plate
            .union(stages)
            .faces("<Z")
            .workplane(origin=(0, offset_y, 0))
            .pushPoints(pts)
            .cboreHole(diameter=thru_dia, cboreDiameter=cbr_dia, cboreDepth=cbr_depth))

def add_internal_bosses(base: Workplane, lid: Workplane,
                        width, length,
                        height=HEIGHT,
                        bottom_thick=BASE_BOTTOM_THICK,
                        base_wall_thick=BASE_WALL_THICK,
                        lid_wall_thick=LID_WALL_THICK,
                        thru_dia=LID_SCREW_HOLE_DIAMETER,
                        boss_dia=BOSS_DIAMETER,
                        nut_dia=NUT_HOLE_DIAMETER,
                        cbr_dia=LID_SCREW_CBR_DIAMETER,
                        cbr_depth=LID_SCREW_CBR_DEPTH) -> tuple[Workplane, Workplane]:
    """
    Add internal bosses in base corners for lid screws
    """
    melt_depth = 1.0
    hx = width / 2.0 - lid_wall_thick - boss_dia / 2.0 + melt_depth
    hy = length / 2.0 - lid_wall_thick - boss_dia / 2.0
    pts = [(hx, hy), (-hx, hy), (hx, -hy + base_wall_thick), (-hx, -hy + base_wall_thick)]
    hole_height = 15.0
    for (x, y) in pts:
        # create a cylinder rising from inside bottom surface
        boss = (cq
                .Workplane("XY")
                .transformed(offset=(x, y, bottom_thick))
                .circle(boss_dia / 2.0)
                .extrude(height - bottom_thick)
                .faces("<Z")
                .workplane(centerOption="CenterOfMass")
                .hole(nut_dia, depth=hole_height)) # hole for screw into boss
        lid = lid.union(boss)

        base = (base
               .faces("<Z")
               .workplane(origin=(0, lid_wall_thick, 0))
               .pushPoints(pts)
               .cboreHole(diameter=thru_dia, cboreDiameter=cbr_dia, cboreDepth=cbr_depth))

    return base, lid

def build_base(width, length,
               height=HEIGHT,
               bottom_thick=BASE_BOTTOM_THICK,
               top_thick=LID_TOP_THICK,
               base_wall_thick=BASE_WALL_THICK,
               lid_wall_thick=LID_WALL_THICK) -> tuple[Workplane, Workplane]:
    """
    Build base: bottom plate + rear wall attached for a covenient mounting of rear connectors
    """
    # bottom plate centered at origin, thickness bottom_thick
    bottom = (cq
              .Workplane("XY")
              .box(width, length, bottom_thick, centered=(True, True, False)))

    rear_length = width - lid_wall_thick
    wall_height = height - bottom_thick - top_thick / 2.0 + bottom_thick / 2.0

    # rear wall
    rear = (cq
            .Workplane("XY")
            .transformed(offset=(0, -length / 2.0 + base_wall_thick / 2.0 + lid_wall_thick, bottom_thick / 2.0))
            .box(rear_length, base_wall_thick, wall_height, centered=(True, True, False)))

    return bottom, rear

def build_lid(width, length,
              height=HEIGHT,
              top_thick=LID_TOP_THICK,
              bottom_thick=BASE_BOTTOM_THICK,
              wall_thick=LID_WALL_THICK) -> tuple[Workplane, Workplane, Workplane, Workplane]:
    """
    Build lid: top plate + 3 side walls (front, left, right)
    """
    # top plate centered at origin but later positioned at top; we'll create centered and translate when assembling
    top = (cq
           .Workplane("XY")
           .transformed(offset=(0, 0, height - top_thick))
           .box(width, length, top_thick, centered=(True, True, False)))

    wall_height = height - bottom_thick / 2.0
    wall_z_offset = bottom_thick / 2.0

    # front wall (positive Y direction)
    front = (cq
             .Workplane("XY")
             .transformed(offset=(0, length / 2.0 - wall_thick / 2.0, wall_z_offset))
             .box(width - wall_thick * 2.0, wall_thick, wall_height, centered=(True, True, False)))

    # left side wall (negative X): create wall and position so it attaches to left edge of top
    left = (cq
            .Workplane("XY")
            .transformed(offset=(-width / 2.0 + wall_thick / 2.0, 0, wall_z_offset))
            .box(wall_thick, length, wall_height, centered=(True, True, False)))

    # right side wall (positive X)
    right = (cq
             .Workplane("XY")
             .transformed(offset=(width / 2.0 - wall_thick / 2.0, 0, wall_z_offset))
             .box(wall_thick, length, wall_height, centered=(True, True, False)))

    return top, front, left, right

def assemble_parts(base_bottom: Workplane, base_rear: Workplane,
                           lid_top: Workplane, lid_front: Workplane,
                           lid_left: Workplane, lid_right: Workplane,
                           fillet_radius: float = FILLET_RADIUS) -> tuple[Workplane, Workplane]:

    base_bottom = (base_bottom
                   .cut(lid_front)
                   .cut(lid_left)
                   .cut(lid_right))
    lid_top = lid_top.cut(base_rear)
    lid_left = lid_left.cut(base_rear)
    lid_right = lid_right.cut(base_rear)

    # Filleting is here for aesthetic and better parts fitting purposes.
    # The following approach allows for some imperfect issue in the rear part of the base,
    # but generally gives good results without making the process yet more complicated.

    # For the sake of simplicity, let's fillet base's parts separately
    base_rear = base_rear.edges(">Z or |Z").fillet(fillet_radius)
    base_bottom = base_bottom.edges("|Z or >Z or <Z").fillet(fillet_radius)

    # Fillet edges of the lid when it's fully assembled
    lid = lid_top.union(lid_front).union(lid_left).union(lid_right)
    lid = lid.edges(">Z or ((>Y or <Y or >X or <X) exc <Z)").fillet(fillet_radius)

    return base_bottom.union(base_rear), lid

# Create rear wall cutouts for preamp
def add_preamp_rear_connectors(base: Workplane,
                               height=HEIGHT,
                               preamp_length=PREAMP_LENGTH,
                               cinch_dia=CINCH_DIAMETER,
                               power_dia=POWER_DIAMETER,
                               ground_dia=GND_DIAMETER,
                               base_wall_thick=BASE_WALL_THICK,
                               lid_wall_thick=LID_WALL_THICK,
                               add_labels=True) -> Workplane:
    # Place input pair centered at X = -30, power at X = 0, output pair at X = +30 (approx)
    cinch_interval = 15.0
    cinch_x_offset = cinch_interval + power_dia
    z_offset = height / 2.0
    input_1_x_offset = -cinch_interval / 2.0 - cinch_x_offset
    output_1_x_offset = -cinch_interval / 2.0 + cinch_x_offset
    base = (base
            .faces("<Y")
            .workplane(offset=-lid_wall_thick, centerOption="CenterOfMass")
            .pushPoints([(input_1_x_offset, z_offset), (cinch_interval / 2.0 - cinch_x_offset, z_offset)])
            .hole(cinch_dia, depth=base_wall_thick)  # Input pair
            .pushPoints([(0, z_offset - ground_dia / 2.0 - cinch_dia)])
            .hole(power_dia, depth=base_wall_thick)  # Power in (center)
            .pushPoints([(0, z_offset + power_dia / 2.0)])
            .hole(ground_dia, depth=base_wall_thick) # Ground terminal
            .pushPoints([(output_1_x_offset, z_offset), (cinch_interval / 2.0 + cinch_x_offset, z_offset)])
            .hole(cinch_dia, depth=base_wall_thick)) # Output pair

    if add_labels:
        text_depth = 0.5
        font_size = 4.0
        text_y_offset = -(preamp_length / 2.0 - lid_wall_thick)
        text_input = (cq
                    .Workplane("XZ")
                    .text("Input", font_size, text_depth)
                    .translate((input_1_x_offset + cinch_dia, text_y_offset, z_offset + cinch_dia * 1.3)))
        text_ground = (cq
                    .Workplane("XZ")
                    .text("Ground", font_size, text_depth)
                    .translate((0, text_y_offset, z_offset + power_dia * 1.3)))
        text_output = (cq
                    .Workplane("XZ")
                    .text("Output", font_size, text_depth)
                    .translate((output_1_x_offset + cinch_dia, text_y_offset, z_offset + cinch_dia * 1.3)))
        text_power = (cq
                    .Workplane("XZ")
                    .text("Power", font_size, text_depth)
                    .translate((power_dia * 1.3, text_y_offset, z_offset - ground_dia / 2.0 - cinch_dia + font_size / 2.0)))
        base = base.union(text_input).union(text_ground).union(text_output).union(text_power)

    return base


# Create rear wall cutouts for PSU
def add_psu_rear_connectors(base: Workplane,
                            height=HEIGHT,
                            psu_length=PSU_LENGTH,
                            psu_width=PSU_WIDTH,
                            power_dia=POWER_DIAMETER,
                            socket_width=POWER_SOCKET_WIDTH,
                            socket_height=POWER_SOCKET_HEIGHT,
                            base_wall_thick=BASE_WALL_THICK,
                            lid_wall_thick=LID_WALL_THICK,
                            add_labels: bool=True) -> Workplane:
    hole_x_offset = 9.0
    z_offset = height / 2.0
    base = (base
            .faces("<Y")
            .workplane(offset=-lid_wall_thick - base_wall_thick, centerOption="CenterOfMass")
            .moveTo(-psu_width / 2.0 + socket_width / 2.0 + hole_x_offset + base_wall_thick, z_offset)
            .rect(socket_width, socket_height)
            .cutBlind(until=base_wall_thick)          # Power socket cutout
            .moveTo(-psu_width / 2.0 + socket_width / 2.0 + hole_x_offset + base_wall_thick, z_offset)
            .rect(socket_width, socket_height + 4.0)
            .cutBlind(until=base_wall_thick-2.0)      # Power socket mounts
            .transformed(offset=(0, 0, lid_wall_thick))
            .pushPoints([(psu_width / 2.0 - power_dia / 2.0 - hole_x_offset - base_wall_thick, z_offset)])
            .hole(power_dia, depth=base_wall_thick))  # Power out

    if add_labels:
        text_depth = 0.5
        font_size = 4.0
        text_y_offset = -(psu_length / 2.0 - lid_wall_thick)
        text_power = (cq
                    .Workplane("XZ")
                    .text("PWR Out", font_size, text_depth)
                    .translate((psu_width / 2.0 - power_dia / 2.0 - hole_x_offset - base_wall_thick, text_y_offset, z_offset + power_dia)))
        base = base.union(text_power)

    return base

# Create bottom and front cutouts for preamp
def add_preamp_control_cutouts(base: Workplane, offset_y,
                               bottom_thick=BASE_BOTTOM_THICK) -> Workplane:
    # Dimensions of the switches' cutouts
    mm_mc_x_len = 12.0
    mm_mc_y_len = 8.0
    cap_res_x_len = 16.0
    cap_res_y_len = 12.0
    # Offsets of the switches' centers from the center of PCB
    cap_res_x_offset = 26.0
    cap_res_y_offset = 49.0
    mm_mc_x_offset = 30.5
    mm_mc_y_offset = 11.0

    return (base
            .faces("<Z")
            .workplane(origin=(0, offset_y, 0))  # Place the origin right in the center of PCB
            .pushPoints([(mm_mc_x_offset, mm_mc_y_offset), (-mm_mc_x_offset, mm_mc_y_offset)])
            .rect(mm_mc_x_len, mm_mc_y_len)
            .cutBlind(until=-bottom_thick)    # Right+left channel MM/MC switch
            .pushPoints([(cap_res_x_offset, -cap_res_y_offset), (-cap_res_x_offset, -cap_res_y_offset)])
            .rect(cap_res_x_len, cap_res_y_len)
            .cutBlind(until=-bottom_thick))   # Right+left channel capacity/resistance switch

if __name__ == "__main__":
    # Build preamp base (bottom + rear wall)
    preamp_bottom, preamp_rear = build_base(PREAMP_WIDTH, PREAMP_LENGTH)
    # Build preamp lid (top + side walls)
    preamp_top, preamp_front, preamp_left, preamp_right = build_lid(PREAMP_WIDTH, PREAMP_LENGTH)
    # Make sure lid fits over base
    preamp_base, preamp_lid = assemble_parts(preamp_bottom, preamp_rear, preamp_top, preamp_front, preamp_left, preamp_right)

    # Add PCB mounting stages and holes
    # Move PCB closer to the front wall to make room for rear connectors
    preamp_pcb_offset_y = (PREAMP_LENGTH - PREAMP_PCB_SPACING_Y) / 2.0 - LID_WALL_THICK - BOSS_DIAMETER - PCB_STAGE_DIAMETER - 3.0
    preamp_base = add_pcb_mounts(preamp_base, PREAMP_PCB_SPACING_X, PREAMP_PCB_SPACING_Y, preamp_pcb_offset_y)

    # Add connector cutouts
    preamp_base = add_preamp_rear_connectors(preamp_base)
    # Add bottom cutouts for controls
    preamp_base = add_preamp_control_cutouts(preamp_base, preamp_pcb_offset_y)
    # Add hole for LED indicator on front wall
    led_diameter = 3.2
    preamp_lid = (preamp_lid
                   .faces(">Y")
                   .workplane(centerOption="CenterOfMass")
                   .pushPoints([(0, 0)])
                   .hole(led_diameter, LID_WALL_THICK))
    # Add internal bosses for screws
    preamp_base, preamp_lid = add_internal_bosses(preamp_base, preamp_lid, PREAMP_WIDTH, PREAMP_LENGTH)

    # Don't cut vents thru, only emboss pattern in the lid side walls
    preamp_vent_length = PREAMP_LENGTH - BOSS_DIAMETER * 2.0 - BASE_WALL_THICK - LID_WALL_THICK * 2.0
    preamp_vent_height = HEIGHT * 0.6
    preamp_lid = preamp_lid.faces("<X").workplane(centerOption="CenterOfMass")
    #preamp_lid = add_vertical_wall_vents(preamp_lid, preamp_vent_length, preamp_vent_height, LID_WALL_THICK / 2.0)
    preamp_lid = add_diagonal_wall_vents(preamp_lid, preamp_vent_length, preamp_vent_height, LID_WALL_THICK / 4.0, PREAMP_WIDTH)
    preamp_lid = preamp_lid.faces(">X").workplane(centerOption="CenterOfMass")
    #preamp_lid = add_vertical_wall_vents(preamp_lid, preamp_vent_length, preamp_vent_height, LID_WALL_THICK / 2.0)
    preamp_lid = add_diagonal_wall_vents(preamp_lid, preamp_vent_length, preamp_vent_height, -LID_WALL_THICK / 4.0, -PREAMP_WIDTH, angle=-45.0)

    print("Preamp case built")

    # Build PSU base (bottom + rear wall)
    psu_bottom, psu_rear = build_base(PSU_WIDTH, PSU_LENGTH)
    # Build PSU lid (top + side walls)
    psu_top, psu_front, psu_left, psu_right = build_lid(PSU_WIDTH, PSU_LENGTH)
    # Make sure lid fits over base
    psu_base, psu_lid = assemble_parts(psu_bottom, psu_rear, psu_top, psu_front, psu_left, psu_right)

    # Add PCB mounting stages and holes
    psu_base = add_pcb_mounts(psu_base, PSU_PCB_SPACING_X, PSU_PCB_SPACING_Y, (PSU_LENGTH - PSU_PCB_SPACING_Y) / 4.0 - BASE_WALL_THICK / 2.0)

    # Add connector cutouts
    psu_base = add_psu_rear_connectors(psu_base)
    # Add internal bosses for screws
    psu_base, psu_lid = add_internal_bosses(psu_base, psu_lid, PSU_WIDTH, PSU_LENGTH)

    # Vents only in the lid side walls
    psu_vent_length = PSU_LENGTH - BOSS_DIAMETER * 2.0 - BASE_WALL_THICK - LID_WALL_THICK * 2.0
    psu_vent_height = HEIGHT * 0.6
    psu_lid = psu_lid.faces("<X").workplane(centerOption="CenterOfMass")
    #psu_lid = add_vertical_wall_vents(psu_lid, psu_vent_length, psu_vent_height, LID_WALL_THICK)
    psu_lid = add_diagonal_wall_vents(psu_lid, psu_vent_length, psu_vent_height, LID_WALL_THICK, PSU_WIDTH)
    psu_lid = psu_lid.faces(">X").workplane(centerOption="CenterOfMass")
    #psu_lid = add_vertical_wall_vents(psu_lid, psu_vent_length, psu_vent_height, LID_WALL_THICK)
    psu_lid = add_diagonal_wall_vents(psu_lid, psu_vent_length, psu_vent_height, -LID_WALL_THICK, -PSU_WIDTH, angle=-45.0)

    print("PSU case built")

    all_models = { "preamp_base": preamp_base, "preamp_lid": preamp_lid, "psu_base": psu_base, "psu_lid": psu_lid }

    # Optional export
    cq_utils.export_models(DO_STL_EXPORT, DO_STEP_EXPORT, **all_models)

    # Show in OCP Viewer
    cq_utils.show_models(**all_models)
