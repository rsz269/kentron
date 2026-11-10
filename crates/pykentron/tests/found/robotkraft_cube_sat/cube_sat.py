# From cq-editor run:
# _p = '/home/manu/manuLinux/code/sources.yguel/3D_modeling/cadquery/cube_sat/cube_sat.py'
# exec(compile(open(_p).read(), _p, "exec"))

import inspect
import math as m
import os
import sys
from importlib import reload

from cadquery import importers


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

import cube_sat__parameters
from cq_utils import *

reload(cube_sat__parameters)
from cube_sat__parameters import *

# Warning view_trough, is not the same piece built for an amovible front cover, the inner cylinder that is used to hollow the object is just extended to cut trough the font face.
view_through = False
export_stl_step = True
SEE_HALF = False
export_tol = 0.001
DEBUG = False


# ==============================
#  Computed geometric dimensions
# ==============================
TR_Y = [87.5, 69.8]
delta = TR_Y[0] - TR_Y[1]
delta_p = [0, 0.15, 0.25]
for i in range(3):
    TR_Y.append(TR_Y[-1] - delta + delta_p[i])
STEP_FILES = [
    current_path + "/cube_brick/Cube_Brick__Part_1.step",
    current_path + "/cube_brick/Cube_Brick__Part_2.step",
    current_path + "/cube_brick/Cube_Brick__Part_3.step",
    current_path + "/cube_brick/Cube_Brick__Part_4.step",
    current_path + "/cube_brick/Cube_Brick__Part_5.step",
    current_path + "/cube_brick/Cube_Brick__Part_6.step",
]
rendered = {i: False for i in range(len(TR_Y))}

####################################
####################################
# ACTUAL CODE
####################################


def cut_box(tr_y: float = TR_Y[0]):
    x_length = 120
    y_width = 4
    z_height = 120
    box = (
        cq.Workplane("XY")
        .box(
            x_length,
            y_width,
            z_height,
            centered=True,
        )
        .translate((0, tr_y, 0))
    )
    return box


def sat(idx: int = 0):
    """
    Remove snap to sat1
    """
    sat_file = STEP_FILES[idx]
    sat = importers.importStep(sat_file)
    if DEBUG and not rendered[idx]:
        show_object(sat, name=f"imported sat {idx + 1}", options={"color": "red"})
        rendered[idx] = True
    box = cut_box(tr_y=TR_Y[idx])
    sat = sat.cut(box)
    return sat


sats_gen = [(lambda i=i: sat(i)) for i in range(len(TR_Y))]


def add_smooth_pins(shape_gen=sats_gen[0], debug: bool = False):
    """Add smooth pins to the given shape."""
    DIR = -1  # -1 = towards -Y (outward for the "<Y" faces), +1 = towards +Y
    HEIGHT = (
        CubeSat.pin_smooth_height
    )  # distance of the disc from the face bounding-box center (along Y)
    RATIO = CubeSat.pin_smooth_disc_ratio  # disc diameter = RATIO * bbox diagonal

    s = shape_gen()
    bot_faces = s.faces("<Y")
    lofts = []

    for f in bot_faces:
        bb = f.BoundingBox()
        center = cq.Vector(bb.center.x, bb.center.y, bb.center.z)

        diameter = RATIO * bb.DiagonalLength
        disc_center = center + cq.Vector(0, DIR * HEIGHT, 0)

        # circle parallel to the XZ plane -> its normal is along Y
        disc_wire = cq.Wire.makeCircle(diameter / 2, disc_center, cq.Vector(0, 1, 0))

        # loft between the face outline and the circle
        solid = cq.Solid.makeLoft([f.outerWire(), disc_wire], ruled=False)
        lofts.append(solid)

        if debug:
            print(f"face center={center.toTuple()}, disc diameter={diameter:.2f}")

    # Display for debug
    if debug:
        show_object(bot_faces, name="selected faces", options={"color": "red"})
        show_object(
            cq.Compound.makeCompound(lofts),
            name="lofts",
            options={"color": "blue", "alpha": 0.3},
        )

    for f in lofts:
        s = s.union(f)
    return s


smooth_pins_sat_gen = [
    (lambda i=i: add_smooth_pins(sats_gen[i])) for i in range(len(TR_Y))
]


def corner_debug_gen(shape_gen=sats_gen[0], debug: bool = False):
    """Just provide a corner for testing the feature."""
    corner_box = (
        cq.Workplane("ZX")
        .box(20, 20, 200, centered=(True, True, True))
        .translate((50, 0, 50))
    )
    s = shape_gen().intersect(corner_box)
    return lambda: s


str_size = "improved"
color_table = ["blue", "green", "red", "yellow", "pink", "cyan"]
full_model_info = {}
for i in range(len(TR_Y)):
    full_model_info[f"sat_part{i + 1}"] = {
        "name": f"sat_part{i + 1}_" + str_size,
        "gen": smooth_pins_sat_gen[i],
        "color": color_table[i % len(color_table)],
        "export": True,
        "display": True,
    }
    if DEBUG:
        full_model_info[f"corner_part{i + 1}"] = {
            "name": f"corner_part{i + 1}_" + str_size,
            "gen": corner_debug_gen(smooth_pins_sat_gen[i]),
            "color": color_table[i % len(color_table)],
            "export": True,
            "display": True,
        }


if export_stl_step:
    import os

    from cadquery import exporters

    ex_path = os.path.join(current_path, "exports", "models")

    sol_pfx = "sat_" + str_size + "_"

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

if SEE_HALF:
    cut_side = 1000
    cut_cube = (
        cq.Workplane("XY")
        .box(cut_side, cut_side, cut_side, centered=True)
        .translate((0, cut_side / 2, 0))
    )


for mod in full_model_info.values():
    if mod["display"]:
        if SEE_HALF:
            show_object(
                mod["gen"]().cut(cut_cube), mod["name"], options={"color": mod["color"]}
            )
        else:
            show_object(mod["gen"](), mod["name"], options={"color": mod["color"]})


if DEBUG:
    # test box cut
    cbs = [cut_box(tr_y=TR_Y[i]) for i in range(len(TR_Y))]
    show_object(cbs[-1], "cut_box", options={"color": "blue", "alpha": 0.5})
    # smooth_pins_sat1()

# Show last piece
last_piece = importers.importStep(STEP_FILES[-1])
show_object(last_piece, name="imported sat 6", options={"color": "red"})


"""
show_object(place_3x_pg9_connector_holes().translate((0, 0, - HexaJoint.height)), "connector_holes_" +
            str_size, options={"color": "red"})

show_object(pg9_connector_hole(), "connector_hole_" +
            str_size, options={"color": "red"})

show_object(clipped_hollow_sphere(50, 50 - 4, 0),
            "clipped_hollow_sphere_" + str_size, options={"color": "yellow"})

show_object(plate_for_cut, "plate_for_cut_" +
            str_size, options={"color": "red"})

cut_for_screws = bottomPlate_fixing_screw_holes_to_base()
show_object(cut_for_screws, "cut_for_screws_" +
            str_size, options={"color": "red"})

show_object(pod_hat, "pod_hat_" + str_size, options={"color": "cyan"})
s_holes = all_screw_holes(
    2 * Base.top_hexagon_height).translate((0, 0, -Base.top_hexagon_height))
show_object(s_holes, "s_holes_" + str_size, options={"color": "red"})
cyl_corner_rad = ShellSection.CylinderCorner.radius
corner_pos_init = (
    ShellSection.Outer.vertex_2_vertex_radius - cyl_corner_rad, 0, 0)
show_object(corner(ShellNav.height).translate(corner_pos_init),
            "corner_" + str_size, options={"color": "red"})
"""
# show_object(corner(10), "corner")
# show_object(rails_orig00(ShellNav.height, 1.6, 3, 2), "rails")
# show_object(rails(ShellNav.height, 1.6, 3, 2)[
#            0].translate((0, 0, ShellNav.height)), "rails2", options={"color": "red"})

if SEE_HALF:
    # show_object(half_base_emptied(), "Main open")
    pass


print("=============")
print(" Cotes utiles")
print("=============")

box_length = full_model_info["sat_part1"]["bbox"]["length"]
box_width = full_model_info["sat_part1"]["bbox"]["width"]
box_height = full_model_info["sat_part1"]["bbox"]["height"]
print(
    f"bbox: x_length={box_length:.1f} mm, y_width={box_width:.1f} mm, z_height={box_height:.1f} mm."
)
