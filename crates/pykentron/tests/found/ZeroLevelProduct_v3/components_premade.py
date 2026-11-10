# ES CORRECTO PERO NO JUNTA TODAS LAS PARTES DEL IHX Y ENTONCES EN EL ONSHAPE SALEN TODOS LOS TUBOS ETC
"""
Pre-made domain components.

Accessed through the same dict interface as build_3D_primitive(), so they
slot into assemble_objects() and build_solid() exactly like any primitive.

Adding a new component
----------------------
1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
2. Add one entry to PREMADE_BUILDERS.
"""

from __future__ import annotations
from typing import Any, cast
import cadquery as cq

from component_premade_reactor_vessel      import create_reactor_vessel
from component_premade_top_plate           import create_top_plate
from component_premade_ihx                 import create_ihx
from component_premade_reactor_core        import create_reactor_core
from component_premade_strongback          import create_strongback
from component_premade_primary_pump        import create_primary_pump
from component_premade_diagrid             import create_diagrid
from component_premade_above_core_structure import create_above_core_structure


def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
    vessel, _ = create_reactor_vessel(
        inner_d            = obj["inner_d"],
        wall_t             = obj["wall_t"],
        straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
        bottom_head_type   = obj.get("bottom_head_type"),
        bottom_head_params = obj.get("bottom_head_params"),
        top_head_type      = obj.get("top_head_type"),
        top_head_params    = obj.get("top_head_params"),
    )
    return vessel


def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
    return create_top_plate(
        plate_outer_d   = obj["outer_d"],
        plate_thickness = obj["thickness"],
        center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
        hole_groups     = obj.get("hole_groups"),
    )


# CORRECT — makeCompound, tubes included and visible
def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
    parts = create_ihx(obj)
    shapes = []
    for s in parts.values():
        if isinstance(s, cq.Workplane):
            shapes.append(s.val())
        else:
            shapes.append(s)
    compound = cq.Compound.makeCompound(shapes)
    return cq.Workplane().newObject([compound])


def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
    return create_reactor_core(
        radius   = obj["radius"],
        height   = obj["height"],
        z_bottom = obj.get("z_bottom", 0.0),
        n_sides  = obj.get("n_sides"),
    )


def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
    return create_strongback(
        total_height            = obj["total_height"],
        flange_radius           = obj["flange_radius"],
        skirt_outer_radius      = obj["skirt_outer_radius"],
        skirt_inner_radius      = obj["skirt_inner_radius"],
        skirt_height            = obj["skirt_height"],
        taper_bottom_z          = obj["taper_bottom_z"],
        bore_radius             = obj["bore_radius"],
        small_hole_radius       = obj["small_hole_radius"],
        small_hole_count        = obj["small_hole_count"],
        small_hole_placement_r  = obj["small_hole_placement_r"],
        z_bottom                = obj.get("z_bottom", 0.0),
        profile_pts             = obj.get("profile_pts"),
    )


def _build_primary_pump(obj: dict[str, Any]) -> cq.Workplane:
    return create_primary_pump(
        barrel_radius  = obj["barrel_radius"],
        barrel_wall_t  = obj["barrel_wall_t"],
        barrel_height  = obj["barrel_height"],
        nozzle_r_pipe  = obj["nozzle_r_pipe"],
        nozzle_wall_t  = obj["nozzle_wall_t"],
        nozzle_L_leg   = obj["nozzle_L_leg"],
        nozzle_R_bend  = obj["nozzle_R_bend"],
        nozzle_arc_deg = obj["nozzle_arc_deg"],
        nozzle_L_inlet = obj["nozzle_L_inlet"],
        nozzle_z       = obj["nozzle_z"],
        flange_width   = obj["flange_width"],
        flange_height  = obj["flange_height"],
        flange_depth   = obj["flange_depth"],
        z_bottom       = obj.get("z_bottom",     0.0),
        flange_z_top   = obj.get("flange_z_top", None),
    )


def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
    """
    Build the hollow diagrid. New kwargs:
      - wall_t_side / wall_t_top / wall_t_bottom (default 0.030 each)
      - open_top, open_bottom (default False)
      - legacy `wall_t` still accepted as an alias for all three walls
    """
    return create_diagrid(
        diameter                = obj["diameter"],
        thickness               = obj["thickness"],
        z_bottom                = obj.get("z_bottom", 0.0),
        wall_t                  = obj.get("wall_t"),
        wall_t_side             = obj.get("wall_t_side",   0.030),
        wall_t_top              = obj.get("wall_t_top",    0.030),
        wall_t_bottom           = obj.get("wall_t_bottom", 0.030),
        open_top                = obj.get("open_top",    False),
        open_bottom             = obj.get("open_bottom", False),
        nozzle_boss_angles_deg  = obj.get("nozzle_boss_angles_deg"),
        nozzle_z_abs            = obj.get("nozzle_z_abs"),
        nozzle_r_bore           = obj.get("nozzle_r_bore",      0.230),
        nozzle_r_boss           = obj.get("nozzle_r_boss",      0.301),
        nozzle_boss_height      = obj.get("nozzle_boss_height", 0.0775),
    )


def _build_above_core_structure(obj: dict[str, Any]) -> cq.Workplane:
    """
    Build the above-core structure (ACS).

    Geometry parameters are forwarded directly to create_above_core_structure().
    Optional `bottom_holes` dict cuts the Ø80/Ø142 hex-pattern through-holes
    + counterbores through the whole structure (including the top cylinder).
    """
    return create_above_core_structure(
        top_cyl_outer_r      = obj["top_cyl_outer_r"],
        top_cyl_height       = obj["top_cyl_height"],
        neck_outer_r         = obj["neck_outer_r"],
        neck_height          = obj["neck_height"],
        collar_outer_r       = obj["collar_outer_r"],
        collar_height        = obj["collar_height"],
        wall_t               = obj["wall_t"],
        cone_bottom_outer_r  = obj["cone_bottom_outer_r"],
        cone_height          = obj["cone_height"],
        bottom_ring_height   = obj["bottom_ring_height"],
        closing_plate_height = obj["closing_plate_height"],
        top_cyl_offset_x     = obj.get("top_cyl_offset_x", 0.0),
        top_cyl_offset_y     = obj.get("top_cyl_offset_y", 0.0),
        flow_hole_groups     = obj.get("flow_hole_groups"),
        bottom_holes         = obj.get("bottom_holes"),
        z_bottom             = obj.get("z_bottom", 0.0),
    )


PREMADE_BUILDERS: dict[str, Any] = {
    "reactor_vessel":        _build_reactor_vessel,
    "reactor_top_plate":     _build_reactor_top_plate,
    "ihx":                   _build_ihx,
    "reactor_core":          _build_reactor_core,
    "strongback":            _build_strongback,
    "primary_pump":          _build_primary_pump,
    "diagrid":                _build_diagrid,
    "above_core_structure":  _build_above_core_structure,
}


def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
    obj_type = obj.get("obj_type", "")
    if obj_type not in PREMADE_BUILDERS:
        raise ValueError(
            f"Unknown premade component {obj_type!r}. "
            f"Available: {sorted(PREMADE_BUILDERS)}"
        )
    return PREMADE_BUILDERS[obj_type](obj)














































# # 19052026 15:20 # ES CORRECTO PERO NO JUNTA TODAS LAS PARTES DEL IHX Y ENTONCES EN EL ONSHAPE SALEN TODOS LOS TUBOS ETC Y NO INCLUYE EL ACS
# """
# Pre-made domain components.

# Accessed through the same dict interface as build_3D_primitive(), so they
# slot into assemble_objects() and build_solid() exactly like any primitive.

# Adding a new component
# ----------------------
# 1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
# 2. Add one entry to PREMADE_BUILDERS.
# """

# from __future__ import annotations
# from typing import Any, cast
# import cadquery as cq

# from component_premade_reactor_vessel import create_reactor_vessel
# from component_premade_top_plate      import create_top_plate
# from component_premade_ihx            import create_ihx
# from component_premade_reactor_core   import create_reactor_core
# from component_premade_strongback     import create_strongback
# from component_premade_primary_pump   import create_primary_pump
# from component_premade_diagrid        import create_diagrid


# def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
#     vessel, _ = create_reactor_vessel(
#         inner_d            = obj["inner_d"],
#         wall_t             = obj["wall_t"],
#         straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
#         bottom_head_type   = obj.get("bottom_head_type"),
#         bottom_head_params = obj.get("bottom_head_params"),
#         top_head_type      = obj.get("top_head_type"),
#         top_head_params    = obj.get("top_head_params"),
#     )
#     return vessel


# def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
#     return create_top_plate(
#         plate_outer_d   = obj["outer_d"],
#         plate_thickness = obj["thickness"],
#         center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
#         hole_groups     = obj.get("hole_groups"),
#     )


# # def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
# #     """
# #     Build the IHX as a single fused solid for STEP export.

# #     create_ihx() returns a dict of 7 sub-components. The tube_bundle is
# #     a cq.Compound of N individual non-touching hollow tube solids — OCCT
# #     cannot fuse non-touching solids into a single solid, so including it
# #     in the fuse loop would keep it as a Compound and produce hundreds of
# #     separate STEP parts in Onshape.

# #     Fix: fuse all structural components (plenums, shell, pipes) into one
# #     solid, and exclude tube_bundle. The bundle_shell (cylindrical envelope)
# #     already represents the tube region visually. Individual tubes are
# #     preserved inside create_ihx() for use in the OpenMC neutronics path.
# #     """
# #     parts = create_ihx(obj)

# #     shapes = []
# #     for key, s in parts.items():
# #         if key == "tube_bundle":
# #             # Non-touching solids — cannot fuse; bundle_shell covers this region.
# #             continue
# #         if isinstance(s, cq.Workplane):
# #             shapes.append(s.val())
# #         else:
# #             shapes.append(s)

# #     if not shapes:
# #         raise ValueError(f"IHX build produced no fuseable shapes for obj_id={obj.get('obj_id')}")

# #     fused = shapes[0]
# #     for s in shapes[1:]:
# #         fused = fused.fuse(s)  # type: ignore

# #     return cq.Workplane().add(fused)


# # CORRECT — makeCompound, tubes included and visible
# def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
#     parts = create_ihx(obj)
#     shapes = []
#     for s in parts.values():
#         if isinstance(s, cq.Workplane):
#             shapes.append(s.val())
#         else:
#             shapes.append(s)
#     compound = cq.Compound.makeCompound(shapes)
#     return cq.Workplane().newObject([compound])







# def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
#     return create_reactor_core(
#         radius   = obj["radius"],
#         height   = obj["height"],
#         z_bottom = obj.get("z_bottom", 0.0),
#         n_sides  = obj.get("n_sides"),
#     )


# def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
#     return create_strongback(
#         total_height            = obj["total_height"],
#         flange_radius           = obj["flange_radius"],
#         skirt_outer_radius      = obj["skirt_outer_radius"],
#         skirt_inner_radius      = obj["skirt_inner_radius"],
#         skirt_height            = obj["skirt_height"],
#         taper_bottom_z          = obj["taper_bottom_z"],
#         bore_radius             = obj["bore_radius"],
#         small_hole_radius       = obj["small_hole_radius"],
#         small_hole_count        = obj["small_hole_count"],
#         small_hole_placement_r  = obj["small_hole_placement_r"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         profile_pts             = obj.get("profile_pts"),
#     )


# def _build_primary_pump(obj: dict[str, Any]) -> cq.Workplane:
#     return create_primary_pump(
#         barrel_radius  = obj["barrel_radius"],
#         barrel_wall_t  = obj["barrel_wall_t"],
#         barrel_height  = obj["barrel_height"],
#         nozzle_r_pipe  = obj["nozzle_r_pipe"],
#         nozzle_wall_t  = obj["nozzle_wall_t"],
#         nozzle_L_leg   = obj["nozzle_L_leg"],
#         nozzle_R_bend  = obj["nozzle_R_bend"],
#         nozzle_arc_deg = obj["nozzle_arc_deg"],
#         nozzle_L_inlet = obj["nozzle_L_inlet"],
#         nozzle_z       = obj["nozzle_z"],
#         flange_width   = obj["flange_width"],
#         flange_height  = obj["flange_height"],
#         flange_depth   = obj["flange_depth"],
#         z_bottom       = obj.get("z_bottom",     0.0),
#         flange_z_top   = obj.get("flange_z_top", None),
#     )


# def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
#     """
#     Build the hollow diagrid. New kwargs:
#       - wall_t_side / wall_t_top / wall_t_bottom (default 0.030 each)
#       - open_top, open_bottom (default False)
#       - legacy `wall_t` still accepted as an alias for all three walls
#     """
#     return create_diagrid(
#         diameter                = obj["diameter"],
#         thickness               = obj["thickness"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         wall_t                  = obj.get("wall_t"),
#         wall_t_side             = obj.get("wall_t_side",   0.030),
#         wall_t_top              = obj.get("wall_t_top",    0.030),
#         wall_t_bottom           = obj.get("wall_t_bottom", 0.030),
#         open_top                = obj.get("open_top",    False),
#         open_bottom             = obj.get("open_bottom", False),
#         nozzle_boss_angles_deg  = obj.get("nozzle_boss_angles_deg"),
#         nozzle_z_abs            = obj.get("nozzle_z_abs"),
#         nozzle_r_bore           = obj.get("nozzle_r_bore",      0.230),
#         nozzle_r_boss           = obj.get("nozzle_r_boss",      0.301),
#         nozzle_boss_height      = obj.get("nozzle_boss_height", 0.0775),
#     )


# PREMADE_BUILDERS: dict[str, Any] = {
#     "reactor_vessel":    _build_reactor_vessel,
#     "reactor_top_plate": _build_reactor_top_plate,
#     "ihx":               _build_ihx,
#     "reactor_core":      _build_reactor_core,
#     "strongback":        _build_strongback,
#     "primary_pump":      _build_primary_pump,
#     "diagrid":           _build_diagrid,
# }


# def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
#     obj_type = obj.get("obj_type", "")
#     if obj_type not in PREMADE_BUILDERS:
#         raise ValueError(
#             f"Unknown premade component {obj_type!r}. "
#             f"Available: {sorted(PREMADE_BUILDERS)}"
#         )
#     return PREMADE_BUILDERS[obj_type](obj)





























































































































# """
# Pre-made domain components.

# Accessed through the same dict interface as build_3D_primitive(), so they
# slot into assemble_objects() and build_solid() exactly like any primitive.

# Adding a new component
# ----------------------
# 1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
# 2. Add one entry to PREMADE_BUILDERS.
# """

# from __future__ import annotations
# from typing import Any, cast
# import cadquery as cq

# from component_premade_reactor_vessel import create_reactor_vessel
# from component_premade_top_plate      import create_top_plate
# from component_premade_ihx            import create_ihx
# from component_premade_reactor_core   import create_reactor_core
# from component_premade_strongback     import create_strongback
# from component_premade_primary_pump   import create_primary_pump
# from component_premade_diagrid        import create_diagrid


# def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
#     vessel, _ = create_reactor_vessel(
#         inner_d            = obj["inner_d"],
#         wall_t             = obj["wall_t"],
#         straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
#         bottom_head_type   = obj.get("bottom_head_type"),
#         bottom_head_params = obj.get("bottom_head_params"),
#         top_head_type      = obj.get("top_head_type"),
#         top_head_params    = obj.get("top_head_params"),
#     )
#     return vessel


# def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
#     return create_top_plate(
#         plate_outer_d   = obj["outer_d"],
#         plate_thickness = obj["thickness"],
#         center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
#         hole_groups     = obj.get("hole_groups"),
#     )


# # def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
# #     parts = create_ihx(obj)
# #     shapes = []
# #     for s in parts.values():
# #         if isinstance(s, cq.Workplane):
# #             shapes.append(s.val())
# #         else:
# #             shapes.append(s)
# #     compound = cq.Compound.makeCompound(shapes)
# #     return cq.Workplane().newObject([compound])


# def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
#     parts = create_ihx(obj)
#     shapes = []
#     for s in parts.values():
#         if isinstance(s, cq.Workplane):
#             shapes.append(s.val())
#         else:
#             shapes.append(s)
#     fused = shapes[0]
#     for s in shapes[1:]:
#         fused = fused.fuse(s)  # type: ignore
#     return cq.Workplane().add(fused)



# def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
#     return create_reactor_core(
#         radius   = obj["radius"],
#         height   = obj["height"],
#         z_bottom = obj.get("z_bottom", 0.0),
#         n_sides  = obj.get("n_sides"),
#     )


# def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
#     return create_strongback(
#         total_height            = obj["total_height"],
#         flange_radius           = obj["flange_radius"],
#         skirt_outer_radius      = obj["skirt_outer_radius"],
#         skirt_inner_radius      = obj["skirt_inner_radius"],
#         skirt_height            = obj["skirt_height"],
#         taper_bottom_z          = obj["taper_bottom_z"],
#         bore_radius             = obj["bore_radius"],
#         small_hole_radius       = obj["small_hole_radius"],
#         small_hole_count        = obj["small_hole_count"],
#         small_hole_placement_r  = obj["small_hole_placement_r"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         profile_pts             = obj.get("profile_pts"),
#     )


# def _build_primary_pump(obj: dict[str, Any]) -> cq.Workplane:
#     return create_primary_pump(
#         barrel_radius  = obj["barrel_radius"],
#         barrel_wall_t  = obj["barrel_wall_t"],
#         barrel_height  = obj["barrel_height"],
#         nozzle_r_pipe  = obj["nozzle_r_pipe"],
#         nozzle_wall_t  = obj["nozzle_wall_t"],
#         nozzle_L_leg   = obj["nozzle_L_leg"],
#         nozzle_R_bend  = obj["nozzle_R_bend"],
#         nozzle_arc_deg = obj["nozzle_arc_deg"],
#         nozzle_L_inlet = obj["nozzle_L_inlet"],
#         nozzle_z       = obj["nozzle_z"],
#         flange_width   = obj["flange_width"],
#         flange_height  = obj["flange_height"],
#         flange_depth   = obj["flange_depth"],
#         z_bottom       = obj.get("z_bottom",     0.0),
#         flange_z_top   = obj.get("flange_z_top", None),
#     )


# def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
#     """
#     Build the hollow diagrid. New kwargs:
#       - wall_t_side / wall_t_top / wall_t_bottom (default 0.030 each)
#       - open_top, open_bottom (default False)
#       - legacy `wall_t` still accepted as an alias for all three walls
#     """
#     return create_diagrid(
#         diameter                = obj["diameter"],
#         thickness               = obj["thickness"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         # wall thicknesses
#         wall_t                  = obj.get("wall_t"),               # legacy
#         wall_t_side             = obj.get("wall_t_side",   0.030),
#         wall_t_top              = obj.get("wall_t_top",    0.030),
#         wall_t_bottom           = obj.get("wall_t_bottom", 0.030),
#         open_top                = obj.get("open_top",    False),
#         open_bottom             = obj.get("open_bottom", False),
#         # boss/bore params (resolver fills these in if you don't)
#         nozzle_boss_angles_deg  = obj.get("nozzle_boss_angles_deg"),
#         nozzle_z_abs            = obj.get("nozzle_z_abs"),
#         nozzle_r_bore           = obj.get("nozzle_r_bore",      0.230),
#         nozzle_r_boss           = obj.get("nozzle_r_boss",      0.301),
#         nozzle_boss_height      = obj.get("nozzle_boss_height", 0.0775),
#     )


# PREMADE_BUILDERS: dict[str, Any] = {
#     "reactor_vessel":    _build_reactor_vessel,
#     "reactor_top_plate": _build_reactor_top_plate,
#     "ihx":               _build_ihx,
#     "reactor_core":      _build_reactor_core,
#     "strongback":        _build_strongback,
#     "primary_pump":      _build_primary_pump,
#     "diagrid":           _build_diagrid,
# }


# def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
#     obj_type = obj.get("obj_type", "")
#     if obj_type not in PREMADE_BUILDERS:
#         raise ValueError(
#             f"Unknown premade component {obj_type!r}. "
#             f"Available: {sorted(PREMADE_BUILDERS)}"
#         )
#     return PREMADE_BUILDERS[obj_type](obj)






































































# """
# Pre-made domain components.

# Accessed through the same dict interface as build_3D_primitive(), so they
# slot into assemble_objects() and build_solid() exactly like any primitive.

# The distinction from components_3D_primitives.py:

#   components_3D_primitives  — pure geometry, no domain knowledge
#                                (cylinder, pipe, box, sphere, …)
#   components_premade        — domain-specific assemblies, built from
#                                primitives + boolean operations
#                                (reactor_vessel, reactor_top_plate, …)

# Adding a new component
# ----------------------
# 1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
# 2. Add one entry to PREMADE_BUILDERS.
# Nothing else in the codebase needs to change.

# Positioning convention
# ----------------------
# Every _build_* function is GEOMETRY ONLY — it builds the solid at or near
# the origin and returns it. Rotation (rotation_angles) and translation
# (center_coords / center_coords_pol) are applied ONCE by build_solid after
# _build_* returns. Never apply positioning inside a _build_* function.
# """

# from __future__ import annotations
# from typing import Any, cast
# import cadquery as cq

# from component_premade_reactor_vessel import create_reactor_vessel
# from component_premade_top_plate      import create_top_plate
# from component_premade_ihx            import create_ihx
# from component_premade_reactor_core   import create_reactor_core
# from component_premade_strongback     import create_strongback
# from component_premade_primary_pump   import create_primary_pump
# from component_premade_diagrid        import create_diagrid


# # ---------------------------------------------------------------------------
# # Individual builders
# # Each takes the full spec dict, extracts geometry params, returns a solid.
# # Positioning (rotation_angles, center_coords) is handled by build_solid —
# # never apply it here.
# # ---------------------------------------------------------------------------

# def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
#     vessel, _ = create_reactor_vessel(
#         inner_d            = obj["inner_d"],
#         wall_t             = obj["wall_t"],
#         straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
#         bottom_head_type   = obj.get("bottom_head_type"),
#         bottom_head_params = obj.get("bottom_head_params"),
#         top_head_type      = obj.get("top_head_type"),
#         top_head_params    = obj.get("top_head_params"),
#     )
#     return vessel


# def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
#     return create_top_plate(
#         plate_outer_d   = obj["outer_d"],
#         plate_thickness = obj["thickness"],
#         center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
#         hole_groups     = obj.get("hole_groups"),
#     )


# def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
#     # create_ihx accepts the full spec dict and reads only the keys it needs.
#     # Extra keys (operation, obj_id, obj_type, center_coords, rotation_angles)
#     # are ignored by create_ihx; positioning is applied by build_solid.
#     parts = create_ihx(obj)
#     shapes = []
#     for s in parts.values():
#         if isinstance(s, cq.Workplane):
#             shapes.append(s.val())
#         else:
#             shapes.append(s)
#     compound = cq.Compound.makeCompound(shapes)
#     return cq.Workplane().newObject([compound])


# def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
#     return create_reactor_core(
#         radius   = obj["radius"],
#         height   = obj["height"],
#         z_bottom = obj.get("z_bottom", 0.0),
#         n_sides  = obj.get("n_sides"),
#     )


# def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
#     return create_strongback(
#         total_height            = obj["total_height"],
#         flange_radius           = obj["flange_radius"],
#         skirt_outer_radius      = obj["skirt_outer_radius"],
#         skirt_inner_radius      = obj["skirt_inner_radius"],
#         skirt_height            = obj["skirt_height"],
#         taper_bottom_z          = obj["taper_bottom_z"],
#         bore_radius             = obj["bore_radius"],
#         small_hole_radius       = obj["small_hole_radius"],
#         small_hole_count        = obj["small_hole_count"],
#         small_hole_placement_r  = obj["small_hole_placement_r"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         profile_pts             = obj.get("profile_pts"),
#     )


# def _build_primary_pump(obj: dict[str, Any]) -> cq.Workplane:
#     # Geometry only — positioning (rotation_angles, center_coords) is applied
#     # by build_solid after this function returns, exactly like every other builder.
#     return create_primary_pump(
#         barrel_radius  = obj["barrel_radius"],
#         barrel_wall_t  = obj["barrel_wall_t"],
#         barrel_height  = obj["barrel_height"],
#         nozzle_r_pipe  = obj["nozzle_r_pipe"],
#         nozzle_wall_t  = obj["nozzle_wall_t"],
#         nozzle_L_leg   = obj["nozzle_L_leg"],
#         nozzle_R_bend  = obj["nozzle_R_bend"],
#         nozzle_arc_deg = obj["nozzle_arc_deg"],
#         nozzle_L_inlet = obj["nozzle_L_inlet"],
#         nozzle_z       = obj["nozzle_z"],
#         flange_width   = obj["flange_width"],
#         flange_height  = obj["flange_height"],
#         flange_depth   = obj["flange_depth"],
#         z_bottom       = obj.get("z_bottom",     0.0),
#         flange_z_top   = obj.get("flange_z_top", None),
#     )


# # def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
# #     return create_diagrid(
# #         diameter  = obj["diameter"],
# #         thickness = obj["thickness"],
# #         z_bottom  = obj.get("z_bottom", 0.0),
# #     )
 
# # def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
# #     return create_diagrid(
# #         diameter           = obj["diameter"],
# #         thickness          = obj["thickness"],
# #         z_bottom           = obj.get("z_bottom", 0.0),
# #         # ── new nozzle boss params ──
# #         pump_angles_deg    = obj.get("pump_angles_deg"),
# #         nozzle_z_abs       = obj.get("nozzle_z_abs"),
# #         nozzle_r_bore      = obj.get("nozzle_r_bore",      0.230),
# #         nozzle_depth       = obj.get("nozzle_depth",       0.300),
# #         nozzle_r_boss      = obj.get("nozzle_r_boss",      0.301),
# #         nozzle_boss_height = obj.get("nozzle_boss_height", 0.080),
# #     )

# def _build_diagrid(obj: dict[str, Any]) -> cq.Workplane:
#     return create_diagrid(
#         diameter                = obj["diameter"],
#         thickness               = obj["thickness"],
#         z_bottom                = obj.get("z_bottom", 0.0),
#         nozzle_boss_angles_deg  = obj.get("nozzle_boss_angles_deg"),  # ← updated
#         nozzle_z_abs            = obj.get("nozzle_z_abs"),
#         nozzle_r_bore           = obj.get("nozzle_r_bore",      0.230),
#        # nozzle_depth            = obj.get("nozzle_depth",       0.300),
#         nozzle_r_boss           = obj.get("nozzle_r_boss",      0.301),
#         nozzle_boss_height      = obj.get("nozzle_boss_height", 0.080),
#     )



# # ---------------------------------------------------------------------------
# # Registry — the only thing that needs editing when adding a new component
# # ---------------------------------------------------------------------------
 
# PREMADE_BUILDERS: dict[str, Any] = {
#     "reactor_vessel":    _build_reactor_vessel,
#     "reactor_top_plate": _build_reactor_top_plate,
#     "ihx":               _build_ihx,
#     "reactor_core":      _build_reactor_core,
#     "strongback":        _build_strongback,
#     "primary_pump":      _build_primary_pump,
#     "diagrid":           _build_diagrid,
# }
 
 
# # ---------------------------------------------------------------------------
# # Public entry point — mirrors build_3D_primitive() signature
# # ---------------------------------------------------------------------------
 
# def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
#     obj_type = obj.get("obj_type", "")
#     if obj_type not in PREMADE_BUILDERS:
#         raise ValueError(
#             f"Unknown premade component {obj_type!r}. "
#             f"Available: {sorted(PREMADE_BUILDERS)}"
#         )
#     return PREMADE_BUILDERS[obj_type](obj)

















































































"""
from __future__ import annotations
from typing import Any, cast
import cadquery as cq

from component_premade_reactor_vessel import create_reactor_vessel
from component_premade_top_plate      import create_top_plate
from component_premade_ihx            import create_ihx
from component_premade_reactor_core   import create_reactor_core
from component_premade_strongback     import create_strongback
from component_premade_primary_pump   import create_primary_pump
from utils import rotate_rpy_about_self_global_axes, move_center_to, convert_polar_to_cartesian


# ---------------------------------------------------------------------------
# Individual builders — each takes the raw dict, returns cq.Workplane
# ---------------------------------------------------------------------------

def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
    vessel, _ = create_reactor_vessel(
        inner_d            = obj["inner_d"],
        wall_t             = obj["wall_t"],
        straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
        bottom_head_type   = obj.get("bottom_head_type"),
        bottom_head_params = obj.get("bottom_head_params"),
        top_head_type      = obj.get("top_head_type"),
        top_head_params    = obj.get("top_head_params"),
    )
    return vessel


def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
    return create_top_plate(
        plate_outer_d   = obj["outer_d"],
        plate_thickness = obj["thickness"],
        center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
        hole_groups     = obj.get("hole_groups"),
    )


def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
    parts = create_ihx(obj)
    shapes = []
    for s in parts.values():
        if isinstance(s, cq.Workplane):
            shapes.append(s.val())
        else:
            shapes.append(s)
    compound = cq.Compound.makeCompound(shapes)
    return cq.Workplane().newObject([compound])


def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
    return create_reactor_core(
        radius   = obj["radius"],
        height   = obj["height"],
        z_bottom = obj.get("z_bottom", 0.0),
        n_sides  = obj.get("n_sides"),
    )


def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
    return create_strongback(
        profile_pts             = obj.get("profile_pts"),
        bore_radius             = obj.get("bore_radius",             0.303),
        small_hole_radius       = obj.get("small_hole_radius",       0.0755),
        small_hole_count        = obj.get("small_hole_count",        6),
        small_hole_placement_r  = obj.get("small_hole_placement_r",  0.900),
        small_hole_z_bottom     = obj.get("small_hole_z_bottom",     0.436),
        small_hole_z_top        = obj.get("small_hole_z_top",        1.242),
        z_bottom                = obj.get("z_bottom",                0.0),
    )


def _build_primary_pump(obj: dict[str, Any]) -> cq.Workplane:
    """
    # Build a primary pump and apply positioning (center_coords / center_coords_pol
    # / rotation_angles) using the same convention as build_solid — so pumps slot
    # into assemble_objects() like any other premade component.
"""
    solid = create_primary_pump(
        barrel_radius  = obj.get("barrel_radius",  1.350 / 2),
        barrel_wall_t  = obj.get("barrel_wall_t",  0.040),
        barrel_height  = obj.get("barrel_height",  12.000),
        nozzle_r_pipe  = obj.get("nozzle_r_pipe",  0.460 / 2),
        nozzle_wall_t  = obj.get("nozzle_wall_t",  0.025),
        nozzle_L_leg   = obj.get("nozzle_L_leg",   0.800),
        nozzle_R_bend  = obj.get("nozzle_R_bend",  0.800),
        nozzle_arc_deg = obj.get("nozzle_arc_deg", 112.5),
        nozzle_z       = obj.get("nozzle_z",       0.450),
        flange_width   = obj.get("flange_width",   0.548),
        flange_height  = obj.get("flange_height",  0.900),
        flange_depth   = obj.get("flange_depth",   0.500),
        flange_z_top   = obj.get("flange_z_top",   11.500),
        z_bottom       = obj.get("z_bottom",       0.0),
    )

    # --- positioning (mirrors build_solid convention) ---
    rotation_angles   = obj.get("rotation_angles",   (0.0, 0.0, 0.0))
    center_coords     = obj.get("center_coords")
    center_coords_pol = obj.get("center_coords_pol")

    roll, pitch, yaw = rotation_angles
    solid = rotate_rpy_about_self_global_axes(solid, roll, pitch, yaw)  # type: ignore

    if center_coords_pol is not None:
        center_coords = convert_polar_to_cartesian(*center_coords_pol)

    if center_coords is not None:
        solid = move_center_to(solid, center_coords)

    return solid


# ---------------------------------------------------------------------------
# Registry — the only thing that needs editing when adding a new component
# ---------------------------------------------------------------------------

PREMADE_BUILDERS: dict[str, Any] = {
    "reactor_vessel":    _build_reactor_vessel,
    "reactor_top_plate": _build_reactor_top_plate,
    "ihx":               _build_ihx,
    "reactor_core":      _build_reactor_core,
    "strongback":        _build_strongback,
    "primary_pump":      _build_primary_pump,
}


# ---------------------------------------------------------------------------
# Public entry point — mirrors build_3D_primitive() signature
# ---------------------------------------------------------------------------

def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
    obj_type = obj.get("obj_type", "")
    if obj_type not in PREMADE_BUILDERS:
        raise ValueError(
            f"Unknown premade component {obj_type!r}. "
            f"Available: {sorted(PREMADE_BUILDERS)}"
        )
    return PREMADE_BUILDERS[obj_type](obj)

"""






























"""
Pre-made domain components.

Accessed through the same dict interface as build_3D_primitive(), so they
slot into assemble_objects() and build_solid() exactly like any primitive.

The distinction from components_3D_primitives.py:

  components_3D_primitives  — pure geometry, no domain knowledge
                               (cylinder, pipe, box, sphere, …)
  components_premade        — domain-specific assemblies, built from
                               primitives + boolean operations
                               (reactor_vessel, reactor_top_plate, …)

Adding a new component
----------------------
1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
2. Add one entry to PREMADE_BUILDERS.
Nothing else in the codebase needs to change.

Usage (identical style to all other build_solid primitives)
-----------------------------------------------------------
>>> RPV = {
...     "operation":          "primitive",
...     "obj_id":             "rpv",
...     "obj_type":           "reactor_vessel",
...     "inner_d":            4.72,
...     "wall_t":             0.04,
...     "straight_h":         5.5,
...     "bottom_head_type":   "ellipsoidal",
...     "bottom_head_params": {"head_depth": 1.0},
... }
>>> TOP_PLATE = {
...     "operation":  "primitive",
...     "obj_id":     "top_plate",
...     "obj_type":   "reactor_top_plate",
...     "outer_d":    4.72 + 2 * 0.04,
...     "thickness":  0.1,
...     "z_bottom":   5.5,
...     "hole_groups": [...],
... }
>>> CORE = {
...     "operation":             "primitive",
...     "obj_id":                "core",
...     "obj_type":              "reactor_core",
...     "radius":                1.8,
...     "height":                3.91,
...     "z_bottom":              1.242,
... }
>>> STRONGBACK = {
...     "operation": "primitive",
...     "obj_id":    "strongback",
...     "obj_type":  "strongback",
...     "z_bottom":  -1.702,
... }
>>> assembly = assemble_objects([RPV, TOP_PLATE, CORE, STRONGBACK])
"""
"""
from __future__ import annotations
from typing import Any, cast
import cadquery as cq

from component_premade_reactor_vessel import create_reactor_vessel
from component_premade_top_plate      import create_top_plate
from component_premade_ihx            import create_ihx
from component_premade_reactor_core   import create_reactor_core
from component_premade_strongback     import create_strongback


# ---------------------------------------------------------------------------
# Individual builders — each takes the raw dict, returns cq.Workplane
# ---------------------------------------------------------------------------

def _build_reactor_vessel(obj: dict[str, Any]) -> cq.Workplane:
    vessel, _ = create_reactor_vessel(
        inner_d            = obj["inner_d"],
        wall_t             = obj["wall_t"],
        straight_h         = cast(float, obj.get("straight_h") or obj.get("height")),
        bottom_head_type   = obj.get("bottom_head_type"),
        bottom_head_params = obj.get("bottom_head_params"),
        top_head_type      = obj.get("top_head_type"),
        top_head_params    = obj.get("top_head_params"),
    )
    return vessel


def _build_reactor_top_plate(obj: dict[str, Any]) -> cq.Workplane:
    return create_top_plate(
        plate_outer_d   = obj["outer_d"],
        plate_thickness = obj["thickness"],
        center_coords   = (0.0, 0.0, obj["z_bottom"] + obj["thickness"] / 2.0),
        hole_groups     = obj.get("hole_groups"),
    )


def _build_ihx(obj: dict[str, Any]) -> cq.Workplane:
    parts = create_ihx(obj)
    shapes = []
    for s in parts.values():
        if isinstance(s, cq.Workplane):
            shapes.append(s.val())
        else:
            shapes.append(s)
    compound = cq.Compound.makeCompound(shapes)
    return cq.Workplane().newObject([compound])


def _build_reactor_core(obj: dict[str, Any]) -> cq.Workplane:
    return create_reactor_core(
        radius   = obj["radius"],
        height   = obj["height"],
        z_bottom = obj.get("z_bottom", 0.0),
        n_sides  = obj.get("n_sides"),
    )


def _build_strongback(obj: dict[str, Any]) -> cq.Workplane:
    return create_strongback(
        profile_pts             = obj.get("profile_pts"),
        bore_radius             = obj.get("bore_radius",             0.303),
        small_hole_radius       = obj.get("small_hole_radius",       0.0755),
        small_hole_count        = obj.get("small_hole_count",        6),
        small_hole_placement_r  = obj.get("small_hole_placement_r",  0.900),
        small_hole_z_bottom     = obj.get("small_hole_z_bottom",     0.436),
        small_hole_z_top        = obj.get("small_hole_z_top",        1.242),
        z_bottom                = obj.get("z_bottom",                0.0),
    )


# ---------------------------------------------------------------------------
# Registry — the only thing that needs editing when adding a new component
# ---------------------------------------------------------------------------

PREMADE_BUILDERS: dict[str, Any] = {
    "reactor_vessel":    _build_reactor_vessel,
    "reactor_top_plate": _build_reactor_top_plate,
    "ihx":               _build_ihx,
    "reactor_core":      _build_reactor_core,
    "strongback":        _build_strongback,
}


# ---------------------------------------------------------------------------
# Public entry point — mirrors build_3D_primitive() signature
# ---------------------------------------------------------------------------

def build_premade_component(obj: dict[str, Any]) -> cq.Workplane:
    obj_type = obj.get("obj_type", "")
    if obj_type not in PREMADE_BUILDERS:
        raise ValueError(
            f"Unknown premade component {obj_type!r}. "
            f"Available: {sorted(PREMADE_BUILDERS)}"
        )
    return PREMADE_BUILDERS[obj_type](obj)
"""


"""
Pre-made domain components.

Accessed through the same dict interface as build_3D_primitive(), so they
slot into assemble_objects() and build_solid() exactly like any primitive.

The distinction from components_3D_primitives.py:

  components_3D_primitives  — pure geometry, no domain knowledge
                               (cylinder, pipe, box, sphere, …)
  components_premade        — domain-specific assemblies, built from
                               primitives + boolean operations
                               (reactor_vessel, reactor_top_plate, …)

Adding a new component
----------------------
1. Write a  _build_<name>(obj: dict) -> cq.Workplane  function below.
2. Add one entry to PREMADE_BUILDERS.
Nothing else in the codebase needs to change.

Usage (identical style to all other build_solid primitives)
-----------------------------------------------------------
>>> RPV = {
...     "operation":          "primitive",
...     "obj_id":             "rpv",
...     "obj_type":           "reactor_vessel",
...     "inner_d":            4.72,
...     "wall_t":             0.04,
...     "straight_h":         5.5,
...     "bottom_head_type":   "ellipsoidal",
...     "bottom_head_params": {"head_depth": 1.0},
... }
>>> TOP_PLATE = {
...     "operation":  "primitive",
...     "obj_id":     "top_plate",
...     "obj_type":   "reactor_top_plate",
...     "outer_d":    4.72 + 2 * 0.04,
...     "thickness":  0.1,
...     "z_bottom":   5.5,
...     "hole_groups": [...],
... }
>>> CORE = {
...     "operation":             "primitive",
...     "obj_id":                "core",
...     "obj_type":              "reactor_core",
...     "radius":                1.8,
...     "height":                3.91,
...     "z_bottom":              1.242,
... }
>>> STRONGBACK = {
...     "operation": "primitive",
...     "obj_id":    "strongback",
...     "obj_type":  "strongback",
...     "z_bottom":  -1.702,
... }
>>> PUMP = {
...     "operation":       "primitive",
...     "obj_id":          "pump_1",
...     "obj_type":        "primary_pump",
...     "barrel_height":   7.0,
...     "z_bottom":        1.0,
...     "center_coords":   (3.369, 0.0, 4.5),
... }
>>> assembly = assemble_objects([RPV, TOP_PLATE, CORE, STRONGBACK, PUMP])
"""