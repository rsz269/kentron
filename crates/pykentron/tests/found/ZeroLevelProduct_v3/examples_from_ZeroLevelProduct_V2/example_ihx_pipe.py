import cadquery as cq
from ocp_vscode import show
from build_3D_solid import build_solid
from assemble import assemble_objects

R_bend  = 1050     # adjust once you have the real value
L_vert  = 8000
L_horiz = 3000
#L_horiz = 6355 - R_bend   # horizontal straight after the bend centre

path_wire = (
    cq.Workplane("XZ")
    .moveTo(0, 0)
    .lineTo(0, L_vert)
    .radiusArc((R_bend, L_vert + R_bend), R_bend)
    .lineTo(R_bend + L_horiz, L_vert + R_bend)
    .wire()
    .val()
)

# path_wire = (
#     cq.Workplane("XZ")
#     .moveTo(0, 0)
#     .lineTo(0, L_vert - 10)       # stop just before the bend
#     .lineTo(0, L_vert)            # short tangent lead-in
#     .radiusArc((R_bend, L_vert + R_bend), -R_bend)
#     .lineTo(R_bend + L_horiz, L_vert + R_bend)
#     .wire()
#     .val()
# )


pipe_436 = {"operation"      : "sweep",
            "profile"        : {"obj_type": "circle", "radius": 436/2},
            "path"           : path_wire,   # type: ignore
            "isFrenet"       : True,
            "wall_thickness" : 20,          # set actual wall thickness once known
            "plane"          : "XY",
            "obj_id"         : "cold_leg_pipe",      # Check if it's the hot or cold leg
            "insert_into"    : "inner_structure",
            #"center_coords" : (0,0,0),
            }

z_pos_inner = +3000

inner = {
     "obj_id" : "inner_structure",
    "operation": "revolve",
    "profile": [(660,0+z_pos_inner),(400,1441+z_pos_inner),(400,5041+z_pos_inner),(0,5041+z_pos_inner),(0,0+z_pos_inner)],
    "plane": "XZ",
    "axis": "Z",
    "angle": 360,
    "wall_thickness": 2,
        
         }


shell = {"obj_id": "outer_shell",
         "operation": "primitive",
         "obj_type": "cylinder_closed_bottom",
         "height": 10000,
         "outer_radius": 663,
         "wall_thickness": 20,
         "center_coords": (0,0,4500),
         "bottom_thickness": 20,
         "bottom_head_type": "hemispherical",
         }


pipe_630 = {"obj_id": "hot_leg_pipe",   # Check if it's the hot or cold leg
            "operation": "primitive",
            "obj_type": "pipe",
            "outer_radius": 630/2,
            "height": 3000,
            "wall_thickness":20,
            "rotation_angles": (0,90,0),
            "center_coords": (1750, 0, 7500),}
#            "insert_into": ["inner_structure", "outer_shell"],}


assembly = assemble_objects([inner, pipe_436, shell, pipe_630])
show(assembly)




# assembly = assemble_objects([pipe_436])
# show(assembly)


#show(pipe_436)


