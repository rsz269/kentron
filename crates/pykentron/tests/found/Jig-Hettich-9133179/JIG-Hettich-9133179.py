"""
Jig design for Hettich 9133179 hinge assembly.

This module generates a CadQuery-based jig for drilling and aligning the Hettich 9133179
hinge components. It creates a base structure with stability ribs and defines methods for
hole creation and positioning.
"""

import cadquery as cq
from cadquery import exporters, Location, Vector

from ocp_vscode import show_object, show

# dimensions in mm
JIG_PLANE_DEPTH = 150
JIG_PLANE_THICKNESS = 10

BOARD_THICKNESS = 28
BOARD_WIDTH = 100
BOARD_LENGTH = 200

HINGE_DIAMETER = 13.5
HINGE_LENGTH = 61.5
HINGE_LENGTH_INNER = 31.7
HINGE_COUNTER_DISTANCE = 46.5

HINGE_DEPTH = 6.5
HINGE_DEPTH_INNER = 18.5

BIT_CUT_LENGTH = 25.4


board_step1 = (
    cq.Workplane("XZ")
    .box(BOARD_LENGTH, BOARD_WIDTH, BOARD_THICKNESS, centered=(True, False, True))
    .translate((0, 0, -BOARD_WIDTH))
)


board_step2 = (
    board_step1.faces(">Z")
    .workplane()
    .slot2D(HINGE_LENGTH, HINGE_DIAMETER, 0)
    .cutBlind(-HINGE_DEPTH)
)


def jig(top_thickness, board):
    """Create a jig for the Hettich 9133179 hinge.

    Returns:
        A CadQuery workplane object representing the complete jig assembly,
        including the base plane, walls, triangular ribs for stability,
        and with the board removed from the structure.
    """

    width = 150
    height = 100
    thickness = BOARD_THICKNESS + 10  # add some extra thickness for wall

    # create the base of the jig
    jig_base = (
        cq.Workplane("XY")
        .rect(width, JIG_PLANE_DEPTH)
        .extrude(-JIG_PLANE_THICKNESS)
        .edges("|Z")
        .fillet(25)
        .faces("<Z")
        .rect(width, thickness)
        .extrude(-height)
        .translate((0, 0, top_thickness))
    )

    # triangular ribs between base plane and wall for stability
    rib_width = JIG_PLANE_DEPTH * 0.75
    rib_depth = 10
    rib_height = rib_width * 0.75
    rib_x_offset = (width / 2) - (rib_width / 2)

    rib = (
        cq.Workplane("YZ")
        .polyline([(-rib_width / 2, 0), (rib_width / 2, 0), (0, -rib_height)])
        .close()
        .extrude(-rib_depth)
    )

    jig_base = jig_base.union(rib.translate((rib_x_offset, 0, top_thickness)))
    jig_base = jig_base.union(rib.translate((-rib_x_offset, 0, top_thickness)))
    jig_base = jig_base - board

    return jig_base


def hinge_top(self, length):
    """Create the top hole for the hinge.

    Args:
        self: The CadQuery workplane object.
        length: The length of the hinge slot.


    Returns:
        The modified workplane with a through-slot cut for the hinge opening,
        sized to HINGE_LENGTH by HINGE_DIAMETER.
    """
    result = self.slot2D(length, HINGE_DIAMETER, 0).cutThruAll()
    return result


cq.Workplane.hinge_top = hinge_top


if __name__ == "__main__":
    TOP_THICKNESS_STEP_1 = BIT_CUT_LENGTH - HINGE_DEPTH + 10  # add some extra thickness
    jig_model = jig(TOP_THICKNESS_STEP_1, board_step1)
    jig_model = jig_model.faces(">Z").workplane().hinge_top(HINGE_LENGTH)

    asm_step1 = (
        cq.Assembly()
        .add(jig_model, name="Jig")
        .add(board_step1, name="Board", color=cq.Color("burlywood1"))
    )

    TOP_THICKNESS_STEP_2 = (
        BIT_CUT_LENGTH - HINGE_DEPTH_INNER + 10
    )  # add some extra thickness
    print(f"Top thickness for step 2: {TOP_THICKNESS_STEP_2} mm")
    jig_model_step2 = jig(TOP_THICKNESS_STEP_2, board_step2)
    jig_model_step2 = jig_model_step2.faces(">Z").workplane().hinge_top(HINGE_LENGTH_INNER)

    asm_step2 = (
        cq.Assembly(loc=Location(0, JIG_PLANE_DEPTH + 20))
        .add(jig_model_step2, name="Jig")
        .add(board_step2, name="Board", color=cq.Color("burlywood1"))
    )

    show_object(asm_step1, name="Step 1")
    show_object(asm_step2, name="Step 2")

    # export the model to a file
    exporters.export(jig_model, "jig-Hettich-9133179.stl")
    exporters.export(jig_model, "jig-Hettich-9133179.step")
    exporters.export(jig_model_step2, "jig-Hettich-9133179_step2.stl")
    exporters.export(jig_model_step2, "jig-Hettich-9133179_step2.step")
