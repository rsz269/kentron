# From cq-editor run:
# _p = '/home/manu/manuLinux/code/sources.yguel/3D_modeling/cadquery/tubes_holder/tubes_holder.py'
# exec(compile(open(_p).read(), _p, "exec"))

import inspect
import math as m
import os
import sys
from importlib import reload


def _script_path():
    if "__file__" in globals():
        return os.path.abspath(__file__)
    filename = inspect.currentframe().f_code.co_filename
    return os.path.abspath(filename)


current_path = os.path.dirname(_script_path())
sys.path.append(current_path)

import cq_utils

reload(cq_utils)
import numpy as np
import scipy.spatial.transform as stf

import tubes_holder__parameters
from cq_utils import *

reload(tubes_holder__parameters)
from tubes_holder__parameters import *

# Warning view_trough, is not the same piece built for an amovible front cover, the inner cylinder that is used to hollow the object is just extended to cut trough the font face.
view_through = False
export_stl_step = True
SEE_HALF = False
export_tol = 0.001


# ==============================
#  Computed geometric dimensions
# ==============================


####################################
####################################
# ACTUAL CODE
####################################

left_foot_translation = (0, 0, 0)
right_foot_translation = (0, 0, 0)


def tube_hole(h: float, tr: tuple = (0, 0, 0)):
    """
    Create a tube hole.
    """
    hole_radius = TubeHolder.tube_hole_radius
    hole = cq.Workplane("XY").circle(hole_radius).extrude(2 * h).translate(tr)
    return hole


def no_stand_holder():
    """
    Create the main holder.
    """
    width = TubeHolder.width
    length = TubeHolder.length
    height = TubeHolder.height
    middle_height = TubeHolder.middle_height
    wall_thickness = TubeHolder.wall_thickness
    base = cq.Workplane("XY").box(
        width, length, wall_thickness, centered=(True, True, False)
    )
    top = base.translate((0, 0, height - wall_thickness))
    middle = base.translate((0, 0, middle_height))
    leg_base = cq.Workplane("XY").box(
        width, wall_thickness, height, centered=(True, True, False)
    )
    left_leg = leg_base.translate((0, -length / 2 + wall_thickness / 2, 0))
    right_leg = leg_base.translate((0, length / 2 - wall_thickness / 2, 0))
    holder = top.union(middle).union(left_leg).union(right_leg)

    # Fillets
    # =======
    fillet_size = wall_thickness / 4.0
    # Add fillets to the edges of the top most face
    holder = holder.edges(">Z").fillet(fillet_size)
    # Add fillets at the junction of the top and middle sections
    ## Left
    l_top_left = (
        width / 2,
        -length / 2 + 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    l_bot_right = (-width / 2, -length / 2 + wall_thickness / 2, wall_thickness / 2)
    left_selector = cq.selectors.BoxSelector(
        l_top_left, l_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(left_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z

    ## Right
    r_top_left = (
        width / 2,
        length / 2 - 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    r_bot_right = (-width / 2, length / 2 - wall_thickness / 2, wall_thickness / 2)
    right_selector = cq.selectors.BoxSelector(
        r_top_left, r_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(right_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z
    return holder


def compute_foot_translations():
    global left_foot_translation, right_foot_translation
    length = TubeHolder.length
    wall_thickness = TubeHolder.wall_thickness
    left_foot_translation = (0, -length / 2 + wall_thickness / 2, 0)
    right_foot_translation = (0, length / 2 - wall_thickness / 2, 0)


def full_stand_holder():
    """
    Create the main holder.
    """
    global left_foot_translation, right_foot_translation
    width = TubeHolder.width
    length = TubeHolder.length
    height = TubeHolder.height
    middle_height = TubeHolder.middle_height
    wall_thickness = TubeHolder.wall_thickness
    base = cq.Workplane("XY").box(
        width, length, wall_thickness, centered=(True, True, False)
    )
    top = base.translate((0, 0, height - wall_thickness))
    middle = base.translate((0, 0, middle_height))
    leg_base = cq.Workplane("XY").box(
        width, wall_thickness, height, centered=(True, True, False)
    )
    stand = cq.Workplane("XY").box(
        2 * width,
        8 * wall_thickness,
        wall_thickness,
        centered=(True, True, False),
    )
    leg_base = leg_base.union(stand)
    compute_foot_translations()
    left_leg = leg_base.translate(left_foot_translation)
    right_leg = leg_base.translate(right_foot_translation)
    holder = top.union(middle).union(left_leg).union(right_leg)

    # Fillets
    # =======
    fillet_size = wall_thickness / 4.0
    # Add fillets to the edges of the top most face
    holder = holder.edges(">Z").fillet(fillet_size)
    # Add fillets at the junction of the top and middle sections
    ## Left
    l_top_left = (
        width / 2,
        -length / 2 + 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    l_bot_right = (-width / 2, -length / 2 + wall_thickness / 2, wall_thickness / 2)
    left_selector = cq.selectors.BoxSelector(
        l_top_left, l_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(left_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z

    ## Right
    r_top_left = (
        width / 2,
        length / 2 - 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    r_bot_right = (-width / 2, length / 2 - wall_thickness / 2, wall_thickness / 2)
    right_selector = cq.selectors.BoxSelector(
        r_top_left, r_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(right_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z
    return holder


def half_stand_holder():
    """
    Create the main holder.
    """
    width = TubeHolder.width
    length = TubeHolder.length
    height = TubeHolder.height
    middle_height = TubeHolder.middle_height
    wall_thickness = TubeHolder.wall_thickness
    base = cq.Workplane("XY").box(
        width, length, wall_thickness, centered=(True, True, False)
    )
    top = base.translate((0, 0, height - wall_thickness))
    middle = base.translate((0, 0, middle_height))
    leg_base = cq.Workplane("XY").box(
        width, wall_thickness, height, centered=(True, True, False)
    )
    stand = (
        cq.Workplane("XY")
        .box(
            1.5 * width,
            8 * wall_thickness,
            wall_thickness,
            centered=(True, True, False),
        )
        .translate((width / 4, 0, 0))
    )
    leg_base = leg_base.union(stand)
    left_leg = leg_base.translate((0, -length / 2 + wall_thickness / 2, 0))
    right_leg = leg_base.translate((0, length / 2 - wall_thickness / 2, 0))
    holder = top.union(middle).union(left_leg).union(right_leg)

    # Fillets
    # =======
    fillet_size = wall_thickness / 4.0
    # Add fillets to the edges of the top most face
    holder = holder.edges(">Z").fillet(fillet_size)
    # Add fillets at the junction of the top and middle sections
    ## Left
    l_top_left = (
        width / 2,
        -length / 2 + 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    l_bot_right = (-width / 2, -length / 2 + wall_thickness / 2, wall_thickness / 2)
    left_selector = cq.selectors.BoxSelector(
        l_top_left, l_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(left_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z

    ## Right
    r_top_left = (
        width / 2,
        length / 2 - 2 * wall_thickness,
        height - wall_thickness / 2,
    )
    r_bot_right = (-width / 2, length / 2 - wall_thickness / 2, wall_thickness / 2)
    right_selector = cq.selectors.BoxSelector(
        r_top_left, r_bot_right, boundingbox=False
    )  # Only center point is ok for edges for bbox inclusion test
    holder = (
        holder.edges(right_selector).edges("#Z").fillet(fillet_size)
    )  # keep only edges orthogonal to Z
    return holder


def tube_holder(holder_gen=no_stand_holder):
    """
    Create the complete tube holder with all tube holes.
    """
    length = TubeHolder.length
    tube_height = TubeHolder.tube_height
    nb_holes = TubeHolder.nb_holes
    hole_diam = TubeHolder.tube_hole_diameter
    spacing = TubeHolder.spacing_between_holes
    width = TubeHolder.width
    height = TubeHolder.height
    wall_thickness = TubeHolder.wall_thickness
    middle_height = TubeHolder.middle_height
    base = holder_gen()
    hole = tube_hole(tube_height, (0, 0, 0))
    for i in range(nb_holes):
        tr = (0, i * (hole_diam + spacing) + 1.5 * hole_diam - length / 2, 0)
        base = base.cut(hole.translate(tr))
    top_selector = cq.selectors.BoxSelector(
        (width / 2, length / 2 + wall_thickness, height + wall_thickness / 2),
        (-width / 2, -length / 2 - wall_thickness, height - wall_thickness / 2),
        boundingbox=True,
    )
    middle_selector = cq.selectors.BoxSelector(
        (width / 2, length / 2 + wall_thickness, middle_height + 2 * wall_thickness),
        (-width / 2, -length / 2 - wall_thickness, middle_height + wall_thickness / 4),
        boundingbox=True,
    )
    base = base.edges(top_selector).edges(">Z").fillet(TubeHolder.fillet_holes)
    base = base.edges(middle_selector).fillet(TubeHolder.fillet_holes)
    return base


def gluable_stand_foot():
    """
    Create a gluable stand foot for the tube holder.
    """
    width = TubeHolder.width
    height = TubeHolder.height
    wall_thickness = TubeHolder.wall_thickness

    leg_base_tol = (
        cq.Workplane("XY")
        .box(
            width + 2 * Tol.fit_e,
            wall_thickness + 2 * Tol.fit_e,
            height,
            centered=(True, True, False),
        )
        .translate((0, 0, wall_thickness))
    )
    stand_up = (
        cq.Workplane("XY")
        .box(
            width + 2 * Tol.fit_e + 2 * wall_thickness,
            3 * wall_thickness + 2 * Tol.fit_e,
            wall_thickness,
            centered=(True, True, False),
        )
        .translate((0, 0, wall_thickness))
    )
    stand = cq.Workplane("XY").box(
        2 * width,
        8 * wall_thickness,
        wall_thickness,
        centered=(True, True, False),
    )

    foot = stand.union(stand_up).cut(leg_base_tol)

    not_filleted_box1 = (
        cq.Workplane("XY")
        .box(
            0.75 * width,
            wall_thickness + 4 * Tol.fit_e,
            height,
            centered=(True, True, False),
        )
        .translate((0, 0, 1.75 * wall_thickness))
    )
    not_filleted_box2 = (
        cq.Workplane("XY")
        .box(
            0.75 * width,
            wall_thickness + 4 * Tol.fit_e,
            height,
            centered=(True, True, False),
        )
        .translate((0, 0, wall_thickness))
    )

    nbox = (not_filleted_box2.cut(not_filleted_box1)).val().BoundingBox()

    edge_selector = -cq.selectors.BoxSelector(
        (nbox.xmin, nbox.ymin, nbox.zmin),
        (nbox.xmax, nbox.ymax, nbox.zmax),
        boundingbox=False,
    )

    # all_edges = foot.edges().vals()
    # box_edges = foot.edges(edge_selector).vals()
    # print(f"All edges: {len(all_edges)}, Box edges: {len(box_edges)}")

    foot = foot.edges(edge_selector).edges("not(<Z)").fillet(TubeHolder.fillet_foot)

    return foot


def full_stand_tube_holder():
    """
    Create the complete tube holder with a full stand.
    """
    return tube_holder(holder_gen=full_stand_holder)


def half_stand_tube_holder():
    """
    Create the complete tube holder with a half stand.
    """
    return tube_holder(holder_gen=half_stand_holder)


str_size = (
    str(TubeHolder.tube_hole_diameter)
    + "x"
    + str(TubeHolder.height)
    + "__nb"
    + str(TubeHolder.nb_holes)
)
full_model_info = {
    "tube_holder": {
        "name": "tube_holder_" + str_size,
        "gen": tube_holder,
        "color": "yellow",
        "export": True,
        "display": True,
    },
    "gluable_stand_foot": {
        "name": "gluable_stand_foot_" + str_size,
        "gen": gluable_stand_foot,
        "color": "orange",
        "export": True,
        "display": False,
    },
    "full_stand_tube_holder": {
        "name": "full_stand_tube_holder_" + str_size,
        "gen": full_stand_tube_holder,
        "color": "blue",
        "export": False,
        "display": False,
    },
    "half_stand_tube_holder": {
        "name": "half_stand_tube_holder_" + str_size,
        "gen": half_stand_tube_holder,
        "color": "green",
        "export": False,
        "display": False,
    },
}


if export_stl_step:
    import os

    from cadquery import exporters

    ex_path = os.path.join(current_path, "exports", "models")

    sol_pfx = "diamHole_x_h__nholes_" + str_size + "_"

    file_names = []
    gen_funcs = []
    ids = []

    for k, mod in full_model_info.items():
        if mod["export"]:
            file_names.append(sol_pfx + k)
            gen_funcs.append(mod["gen"])
            ids.append(k)

    models = []
    for i in range(0, len(file_names)):
        d = {"file_name": file_names[i], "gen_f": gen_funcs[i], "id": ids[i]}
        models.append(d)

    for item in models:
        f_name = item["file_name"]
        gen_f = item["gen_f"]

        shape = gen_f()

        # STL
        stl_full_file_name = os.path.join(ex_path, f_name + ".stl")
        exporters.export(
            shape,
            stl_full_file_name,
            exportType=exporters.ExportTypes.STL,
            tolerance=export_tol,
        )

        # STEP
        step_full_file_name = os.path.join(ex_path, f_name + ".step")
        exporters.export(
            shape,
            step_full_file_name,
            exportType=exporters.ExportTypes.STEP,
            tolerance=export_tol,
        )

        # Compute model data
        bbox = shape.val().BoundingBox()

        full_model_info[item["id"]]["bbox"] = {
            "bbox": bbox,
            "length": bbox.xlen,
            "width": bbox.ylen,
            "height": bbox.zlen,
        }


# show_object(res)
proto_base = False
if proto_base:
    for mod in full_model_info.values():
        if "base_" != mod["name"][0 : len("base_")]:
            mod["display"] = False

cut_cube = None

if SEE_HALF:
    width = TubeHolder.width
    length = TubeHolder.length
    height = TubeHolder.height
    cut_cube = cq.Workplane("XY").box(
        -2 * width, 4 * length, 4 * height, centered=(False, True, True)
    )


for mod in full_model_info.values():
    if mod["display"]:
        if SEE_HALF:
            show_object(
                mod["gen"]().cut(cut_cube), mod["name"], options={"color": mod["color"]}
            )
        else:
            show_object(mod["gen"](), mod["name"], options={"color": mod["color"]})

"""
# test middle selector
width = TubeHolder.width
length = TubeHolder.length
wall_thickness = TubeHolder.wall_thickness
middle_height = TubeHolder.middle_height
top_point = (width / 2, length / 2 + wall_thickness, middle_height + 2 * wall_thickness)
bot_point = (
    -width / 2,
    -length / 2 - wall_thickness,
    middle_height + wall_thickness / 4,
)
box = (
    cq.Workplane("XY")
    .box(
        top_point[0] - bot_point[0],
        top_point[1] - bot_point[1],
        top_point[2] - bot_point[2],
        centered=True,
    )
    .translate(
        (
            (top_point[0] + bot_point[0]) / 2,
            (top_point[1] + bot_point[1]) / 2,
            (top_point[2] + bot_point[2]) / 2,
        )
    )
)

show_object(box, "middle_selector_box", options={"color": "blue"})
"""


def display_feet(cut_shape=None):
    compute_foot_translations()
    if cut_shape is not None:
        left_foot = gluable_stand_foot().translate(left_foot_translation).cut(cut_shape)
        right_foot = (
            gluable_stand_foot().translate(right_foot_translation).cut(cut_shape)
        )
    else:
        wall_thickness = TubeHolder.wall_thickness
        left_foot = (
            gluable_stand_foot()
            .translate(left_foot_translation)
            .translate((0, 0, -wall_thickness))
        )
        right_foot = (
            gluable_stand_foot()
            .translate(right_foot_translation)
            .translate((0, 0, -wall_thickness))
        )
    show_object(left_foot, "left_foot", options={"color": "green"})
    show_object(right_foot, "right_foot", options={"color": "green"})


if SEE_HALF:
    display_feet(cut_shape=cut_cube)
else:
    display_feet()


print("=============")
print(" Cotes utiles")
print("=============")

print(f"tube_diameter: {TubeHolder.tube_diameter:.1f} mm.")
print(f"tube_hole_diameter: {TubeHolder.tube_hole_diameter:.1f} mm.")
print(f"nb_holes: {TubeHolder.nb_holes}")
print(f"width: {TubeHolder.width:.1f} mm.")
print(f"length: {TubeHolder.length:.1f} mm.")
print(f"height: {TubeHolder.height:.1f} mm.")
box_length = full_model_info["tube_holder"]["bbox"]["length"]
box_width = full_model_info["tube_holder"]["bbox"]["width"]
box_height = full_model_info["tube_holder"]["bbox"]["height"]
print(
    f"bbox: length={box_length:.1f} mm, width={box_width:.1f} mm, height={box_height:.1f} mm."
)
