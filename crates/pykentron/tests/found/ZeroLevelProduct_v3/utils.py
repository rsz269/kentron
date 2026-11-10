import cadquery as cq
from typing import Callable, Tuple, Literal, Union, cast
import numpy as np
import logging
import math


logger = logging.getLogger(__name__)

PlaneName = Literal["XY", "XZ", "YZ"]
AxisName = Literal["X", "Y", "Z"]

def export_step(shape, path: str):
    """Export a CadQuery solid or assembly to STEP."""
    if isinstance(shape, cq.Assembly):
        # cq.exporters.export(shape.toCompound(), path, exportType="STEP")
        shape.save(path, exportType="STEP")          # ← preserves names
    else:
        cq.exporters.export(shape, path, exportType="STEP")

def export_stl(shape, path: str, tolerance: float = 0.01):
    """Export a CadQuery solid or assembly to STL."""
    if isinstance(shape, cq.Assembly):
        cq.exporters.export(shape.toCompound(), path, exportType="STL", tolerance=tolerance)
    else:
        cq.exporters.export(shape, path, exportType="STL", tolerance=tolerance)




# -----------------------------------------------------------------------------------------------------------------------------------
# ------------------- FUNCTIONS TO POSITION OBJECTS -----------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------
def convert_polar_to_cartesian(r: float, theta: float, z: float) -> Tuple[float, float, float]:
    """
    Convert cylindrical polar coordinates to Cartesian coordinates.
    
    Args:
        r: radial distance from Z-axis
        theta: azimuthal angle in radians (0 = +X axis)
        z: height
    
    Returns:
        (x, y, z) in Cartesian coordinates
    """
    x = r * math.cos(theta)
    y = r * math.sin(theta)
    return (x, y, z)


# Function used to rotate the objects ---- rotation angles are in DEG ----- X (roll), Y (pitch), Z (yaw) through center
def rotate_rpy_about_self_global_axes(res, roll=0, pitch=0, yaw=0): 
    c = res.val().Center()
    if roll:  res = res.rotate((c.x,c.y,c.z), (c.x+1,c.y,c.z), roll)   # global X through center
    if pitch: res = res.rotate((c.x,c.y,c.z), (c.x,c.y+1,c.z), pitch)  # global Y through center
    if yaw:   res = res.rotate((c.x,c.y,c.z), (c.x,c.y,c.z+1), yaw)    # global Z through center
    return res
    

def move_center_to(res, target_xyz):
    c = res.val().Center()
    dx = target_xyz[0] - c.x
    dy = target_xyz[1] - c.y
    dz = target_xyz[2] - c.z
    return res.translate((dx, dy, dz))





# -----------------------------------------------------------------------------------------------------------------------------------
# -------------------- OPERATIONS ON SINGLE 2D PROFILE TO CREATE A 3D OBJECT BY EXTRUSION/ REVOLUTION/ SWEEP ------------------------
# -----------------------------------------------------------------------------------------------------------------------------------
def extrude_profile(profile_wp: cq.Workplane, height: float, both: bool = False) -> cq.Workplane:
    """Extrude a closed 2D profile. Set both=True for symmetric extrusion."""
    if height == 0: raise ValueError("height must not be 0")
    return profile_wp.extrude(height, both=both)



def revolve_profile(
    profile: cq.Workplane,
    angle: float = 360.0,
    axis: AxisName = "Z",
    axis_point: Tuple[float, float, float] = (0, 0, 0),
) -> cq.Workplane:
    # plane_name = getattr(profile, '_plane_name', 'XY')
    # plane = cq.Workplane(plane_name).plane
    plane = profile.plane

    local_point = cast(cq.Vector, plane.toLocalCoords(cq.Vector(*axis_point)))

    _DIRS = {
        "X": cq.Vector(1, 0, 0),
        "Y": cq.Vector(0, 1, 0),
        "Z": cq.Vector(0, 0, 1),
    }
    if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
    local_dir = cast(cq.Vector, plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin))  # type: ignore[operator]
    # .x and .y here are local plane axes, not global X/Y
    # toLocalCoords handles the mapping for any plane (XY, XZ, YZ)
    
    # p0 = (local_point.x, local_point.y)
    # p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

    # Use whichever two local coordinates carry the axis direction
    if abs(local_dir.z) > abs(local_dir.x) and abs(local_dir.z) > abs(local_dir.y):
        # axis direction is mostly in local Z — use x and z instead
        p0 = (local_point.x, local_point.z)
        p1 = (local_point.x + local_dir.x, local_point.z + local_dir.z)
    else:
        p0 = (local_point.x, local_point.y)
        p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)



    print(f"DEBUG: axis={axis}")
    print(f"  local_dir={local_dir.toTuple()}")
    print(f"  p0={p0}, p1={p1}")


    return profile.revolve(angle, p0, p1)



def sweep_profile(
    profile_wp: cq.Workplane,
    path: Union[
        cq.Wire,
        Callable[[float], Tuple[float, float, float]],
        Tuple[Tuple[float,float,float], ...]
    ],
    num_path_points: int = 50,
    isFrenet: bool = False,
) -> cq.Workplane:
    """
    Sweep a closed 2D profile along a path.

    Path can be provided as:
      1. cq.Wire  — pre-built wire using Wire.make* or Wire.assembleEdges([Edge.make*])
                    e.g. Wire.makeHelix, Wire.assembleEdges([Edge.makeSpline(...)])
                    No origin constraint — user controls positioning.
      2. Callable — analytical function f(t) -> (x, y, z), t in [0, 1]
                    Must start at (0, 0, 0). Approximated as a spline.
      3. Tuple    — ((x0,y0,z0), (x1,y1,z1), ...) two or more points connected
                    by straight line segments. Must start at (0, 0, 0).
    """

    # ----------------------------------------------------------------
    # 1. Build path wire — and track first_edge when we know it
    # ----------------------------------------------------------------
    first_edge = None  # will be set directly when we build the wire ourselves

    if isinstance(path, cq.Wire):
        path_wire = path
        # first_edge resolved below via connectivity search

    elif callable(path):
        path_start = path(0)
        if not np.allclose(path_start, [0, 0, 0], atol=1e-6):
            raise ValueError(f"Analytical path must start at (0,0,0), got {path_start}")
        t_vals     = np.linspace(0, 1, num_path_points)
        path_pts   = [cq.Vector(*path(t)) for t in t_vals]
        path_edge  = cq.Edge.makeSpline(path_pts)
        path_wire  = cq.Wire.assembleEdges([path_edge])
        first_edge = path_edge   # single edge, known

    elif (
        isinstance(path, tuple)
        and len(path) >= 2
        and all(isinstance(p, (tuple, list)) for p in path)
    ):
        if not np.allclose(path[0], [0, 0, 0], atol=1e-6):
            raise ValueError(f"Straight line path must start at (0,0,0), got {path[0]}")
        pts        = [cq.Vector(*p) for p in path]
        edges      = [cq.Edge.makeLine(pts[i], pts[i+1]) for i in range(len(pts) - 1)]
        path_wire  = cq.Wire.assembleEdges(edges)
        first_edge = edges[0]    # known — first edge we built

    else:
        raise TypeError(
            "path must be one of:\n"
            "  - cq.Wire via Wire.make* or Wire.assembleEdges([Edge.make*])\n"
            "  - Callable f(t) -> (x,y,z), t in [0,1], must start at (0,0,0)\n"
            "  - Tuple of 2+ points ((x0,y0,z0), (x1,y1,z1), ...), must start at (0,0,0)"
        )

    # ----------------------------------------------------------------
    # 2. Resolve first_edge for cq.Wire case via connectivity search
    # ----------------------------------------------------------------
    if first_edge is None:
        all_edges = path_wire.Edges()
        if len(all_edges) == 1:
            first_edge = all_edges[0]
        elif path_wire.IsClosed():      # type: ignore[attr-defined]
            first_edge = all_edges[0]   # closed: no unique start, any edge is fine
        else:
            # open wire: the true first edge is the one whose startPoint
            # is not the endPoint of any other edge
            all_end_pts = [e.endPoint() for e in all_edges]     # type: ignore[attr-defined]
            first_edge = next(
                e for e in all_edges
                if not any(
                    np.allclose(e.startPoint().toTuple(), ep.toTuple(), atol=1e-6)
                    for ep in all_end_pts
                )
            )

    # ----------------------------------------------------------------
    # 3. Reposition profile at path start, normal = path tangent
    # ----------------------------------------------------------------
    sketch = profile_wp.val()
    if not isinstance(sketch, cq.Sketch):
        raise TypeError("profile_wp must have a cq.Sketch on the stack")

    start_point = cq.Vector(first_edge.startPoint())
    tangent     = cq.Vector(first_edge.tangentAt(0)).normalized()

    candidates  = [cq.Vector(0, 0, 1), cq.Vector(0, 1, 0), cq.Vector(1, 0, 0)]
    ref         = min(candidates, key=lambda v: abs(tangent.dot(v)))
    x_dir       = (ref - tangent * tangent.dot(ref)).normalized()
    plane       = cq.Plane(origin=start_point, normal=tangent, xDir=x_dir)
    wp          = cq.Workplane(plane).placeSketch(sketch.clean())

    # ----------------------------------------------------------------
    # 4. Sweep
    # ----------------------------------------------------------------
    return wp.sweep(path_wire, isFrenet=isFrenet, makeSolid=True)















# def sweep_profile(
#     profile_wp: cq.Workplane,
#     path: Union[
#         cq.Wire,
#         Callable[[float], Tuple[float, float, float]],
#         Tuple[Tuple[float,float,float], Tuple[float,float,float]]
#     ],
#     num_path_points: int = 50,
#     isFrenet: bool = False,
# ) -> cq.Workplane:
#     """
#     Sweep a closed 2D profile along a path.

#     Path can be provided as:
#       1. cq.Wire  — pre-built wire using Wire.make* or Wire.assembleEdges([Edge.make*])
#                     e.g. Wire.makeHelix, Wire.assembleEdges([Edge.makeSpline(...)])
#                     No origin constraint — user controls positioning.
#       2. Callable — analytical function f(t) -> (x, y, z), t in [0, 1]
#                     Must start at (0, 0, 0).
#       3. Tuple    — ((x0,y0,z0), (x1,y1,z1)) straight line shorthand.
#                     Start must be (0, 0, 0).
#     """

#     # ----------------------------------------------------------------
#     # 1. Build path wire
#     # ----------------------------------------------------------------
#     if isinstance(path, cq.Wire):
#         path_wire = path

#     elif callable(path):
#         path_start = path(0)
#         if not np.allclose(path_start, [0, 0, 0], atol=1e-6):
#             raise ValueError(f"Analytical path must start at (0,0,0), got {path_start}")
#         t_vals    = np.linspace(0, 1, num_path_points)
#         path_pts  = [cq.Vector(*path(t)) for t in t_vals]
#         path_edge = cq.Edge.makeSpline(path_pts)
#         path_wire = cq.Wire.assembleEdges([path_edge])


    

#     # This version is version 2
#     elif (
#     isinstance(path, tuple)
#     and len(path) >= 2
#     and all(isinstance(p, (tuple, list)) for p in path)
#     ):
#         if not np.allclose(path[0], [0, 0, 0], atol=1e-6):
#             raise ValueError(f"Straight line path must start at (0,0,0), got {path[0]}")
#         pts = [cq.Vector(*p) for p in path]
#         path_wire = cq.Wire.assembleEdges([
#             cq.Edge.makeLine(pts[i], pts[i+1])
#             for i in range(len(pts) - 1)
#         ])
    
#     # # This version is version 1
#     # elif (
#     #     isinstance(path, tuple) and len(path) == 2
#     #     and isinstance(path[0], (tuple, list))
#     #     and isinstance(path[1], (tuple, list))
#     # ):
#     #     start, end = path
#     #     if not np.allclose(start, [0, 0, 0], atol=1e-6):
#     #         raise ValueError(f"Straight line path must start at (0,0,0), got {start}")
#     #     path_edge = cq.Edge.makeLine(cq.Vector(*start), cq.Vector(*end))
#     #     path_wire = cq.Wire.assembleEdges([path_edge])

#     else:
#         raise TypeError(
#             "path must be one of:\n"
#             "  - cq.Wire via Wire.make* or Wire.assembleEdges([Edge.make*])\n"
#             "  - Callable f(t) -> (x,y,z), t in [0,1], must start at (0,0,0)\n"
#             "  - Tuple ((x0,y0,z0), (x1,y1,z1)), start must be (0,0,0)"
#         )

#     # ----------------------------------------------------------------
#     # 2. Reposition profile at path start, normal = path tangent
#     # ----------------------------------------------------------------
#     sketch = profile_wp.val()
#     if not isinstance(sketch, cq.Sketch):
#         raise TypeError("profile_wp must have a cq.Sketch on the stack")

#     #start_point = cq.Vector(path_wire.Edges()[0].startPoint())
#     #tangent     = cq.Vector(path_wire.Edges()[0].tangentAt(0)).normalized()
#     start_point = cq.Vector(path_wire.startPoint())
#     tangent     = cq.Vector(path_wire.tangentAt(0)).normalized()
#     world_z = cq.Vector(0, 0, 1)
#     ref     = cq.Vector(1, 0, 0) if abs(tangent.dot(world_z)) > 0.9 else world_z
#     x_dir   = (ref - tangent * tangent.dot(ref)).normalized()
#     plane   = cq.Plane(origin=start_point, normal=tangent, xDir=x_dir)
#     wp      = cq.Workplane(plane).placeSketch(sketch.clean())

#     # ----------------------------------------------------------------
#     # 3. Sweep
#     # ----------------------------------------------------------------
#     return wp.sweep(path_wire, isFrenet=isFrenet, makeSolid=True)



























# -----------------------------------------------------------------------------------------------------------------------------------
# -------------------- OPERATIONS BETWEEN TWO OBJECTS -------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------------------
# def insert_into(base: cq.Workplane, insert: cq.Workplane) -> cq.Workplane:
#     return base.cut(insert).union(insert)


# The following is the old insert_into which did not handle well hollow shapes that were not primitives
# def insert_into(base: cq.Workplane, insert: cq.Workplane) -> cq.Workplane:
#     insert_def = getattr(insert, "_def", {})
#     cutter = getattr(insert, "_outer", insert)  # use outer solid if available, else insert as-is

#     obj_type = insert_def.get("obj_type", "")

#     if obj_type in ("pipe", "cylinder_closed_bottom"):
#         bb = insert.val().BoundingBox()                  # type: ignore[attr-defined]
#         cx = (bb.xmin + bb.xmax) / 2
#         cy = (bb.ymin + bb.ymax) / 2
#         cz = (bb.zmin + bb.zmax) / 2
#         cutter = cq.Workplane("XY").cylinder(insert_def["height"], insert_def["outer_radius"])
#         cutter = move_center_to(cutter, (cx, cy, cz))
#     else:
#         cutter = insert

#     return base.cut(cutter).union(insert)

def insert_into(base: cq.Workplane, insert: cq.Workplane) -> cq.Workplane:
    cutter = getattr(insert, "_outer", insert)  # use outer solid if available, else insert as-is
    return base.cut(cutter).union(insert)














# ESTA FUNCTIONA!!!!!!!!!
# def insert_into(base: cq.Workplane, insert: cq.Workplane) -> cq.Workplane:
#     """
#     Insert an object into a base solid, carving out its occupied volume first.

#     For hollow objects (pipe, cylinder_closed_bottom), cuts the outer solid volume
#     from the base to remove displaced material, then unions the hollow insert back in
#     so the bore remains empty. For solid objects, the insert itself is used as the
#     cutter, which is equivalent to a union.

#     The insert type is read automatically from insert._def, which is attached by
#     build_solid at construction time — no need to pass the definition separately.

#     Args:
#         base:   The solid being cut into (e.g. a wall, block, or vessel).
#         insert: The object being inserted, built with build_solid.

#     Returns:
#         A single cq.Workplane representing the combined solid.

#     Example:
#         >>> pipe, _ = build_solid("primitive", {"obj_type": "pipe", "height": 20, "outer_radius": 3, "inner_radius": 2})
#         >>> wall, _ = build_solid("primitive", {"obj_type": "box", "length": 10, "width": 10, "height": 2})
#         >>> result = insert_into(wall, pipe)
#     """
#     insert_def = getattr(insert, "_def", {})
#     obj_type = insert_def.get("obj_type", "")

#     if obj_type == "pipe":
#         cutter = cq.Workplane("XY").cylinder(insert_def["height"], insert_def["outer_radius"])
#     elif obj_type == "cylinder_closed_bottom":
#         cutter = cq.Workplane("XY").cylinder(insert_def["height"], insert_def["outer_radius"])
#     else:
#         cutter = insert

#     return base.cut(cutter).union(insert)


# def insert_into(base: cq.Workplane, insert: cq.Workplane, insert_def: dict) -> cq.Workplane:
#     """
#     Insert an object into a base solid, carving out its occupied volume first.
#     For hollow objects (e.g. pipes), carves the outer solid, then unions the hollow insert back.
#     For solid objects, the insert itself is the cutter.
#     """
#     if insert_def["obj_type"] == "pipe":
#         cutter = cq.Workplane("XY").cylinder(insert_def["height"], insert_def["outer_radius"])
#     else:
#         cutter = insert
#     return base.cut(cutter).union(insert)













# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     """Revolve a closed 2D profile around a global axis"""
#     wire = profile_wp.val()
#     if not isinstance(wire, cq.Wire):
#         raise TypeError("Expected a closed Wire")
    
#     face = cq.Face.makeFromWires(wire)
#     p0 = cq.Vector(*axis_point)
    
#     if axis == "X":   d = cq.Vector(1, 0, 0)
#     elif axis == "Y": d = cq.Vector(0, 1, 0)
#     elif axis == "Z": d = cq.Vector(0, 0, 1)
#     else: raise ValueError("axis must be 'X', 'Y', or 'Z'")
    
#     p1 = p0 + d
#     logger.debug(f"Revolving profile {angle}° around {axis}-axis")
    
#     solid = cq.Solid.revolve(face.outerWire(), face.innerWires(), angle, p0, p1)
#     return cq.Workplane(obj=solid)


# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_offset: Tuple[float, float] = (0, 0),
# ) -> cq.Workplane:
#     """
#     Revolve a closed 2D profile around a global axis.

#     axis_offset shifts the axis in the two directions perpendicular to it,
#     e.g. for axis="Z", offset=(x, y) moves the Z-axis to (x, y, *).
#     Useful for off-center revolves without needing two explicit points.
#     """
#     if angle == 0: raise ValueError("angle must not be 0")

#     _AXES = {
#         "X": (cq.Vector(1, 0, 0), (1, 2)),  # offset in Y, Z
#         "Y": (cq.Vector(0, 1, 0), (0, 2)),  # offset in X, Z
#         "Z": (cq.Vector(0, 0, 1), (0, 1)),  # offset in X, Y
#     }
#     if axis not in _AXES: raise ValueError("axis must be 'X', 'Y', or 'Z'")

#     direction, (i, j) = _AXES[axis]
#     p0 = cq.Vector(0, 0, 0)
#     p0_coords = [0.0, 0.0, 0.0]
#     p0_coords[i] = axis_offset[0]
#     p0_coords[j] = axis_offset[1]
#     p0 = cq.Vector(*p0_coords)
#     p1 = p0 + direction

#     logger.debug(f"Revolving profile {angle}° around {axis}-axis at offset {axis_offset}")
#     return profile_wp.revolve(angle, p0, p1)


# -----------------------------------------------------------------
# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     """
#     Revolve a closed 2D profile around a global axis.
#     Accepts Wire (polyline), Face, Sketch (built-in 2D sketch), and Compound profiles.
#     """
#     if angle == 0: raise ValueError("angle must not be 0")

#     _DIRS = {"X": cq.Vector(1,0,0), "Y": cq.Vector(0,1,0), "Z": cq.Vector(0,0,1)}
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")

#     val = profile_wp.val()

#     if isinstance(val, cq.Wire):
#         face = cq.Face.makeFromWires(val)
#     elif isinstance(val, cq.Face):
#         face = val
#     elif isinstance(val, cq.Sketch):
#         compound = val._faces
#         wire = cast(cq.Wire, cq.Wire.combine(compound.Wires())[0])
#         face = cq.Face.makeFromWires(wire)
#     elif isinstance(val, cq.Compound):
#         faces = profile_wp.faces().vals()
#         if len(faces) != 1: raise ValueError(f"Expected 1 face, got {len(faces)}")
#         face = faces[0]
#     else:
#         raise TypeError(f"Unsupported stack type: {type(val).__name__}. Expected Wire, Face, Sketch, or Compound.")

#     p0 = cq.Vector(*axis_point)
#     p1 = p0 + _DIRS[axis]

#     solid = cq.Solid.revolve(face.outerWire(), face.innerWires(), angle, p0, p1)  # type: ignore
#     return cq.Workplane(obj=solid)


# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     """
#     Revolve a closed 2D profile around a global axis.
#     Works with Wire-based profiles (create_wire_workplane_straight_connections).
#     axis_point: a point the axis passes through, in global coordinates.
#     """
#     wire = profile_wp.val()
#     if not isinstance(wire, cq.Wire):
#         raise TypeError("Expected a closed Wire on the stack (did you call .close()?)")

#     # Convert global axis_point to local workplane coordinates
#     plane = profile_wp.plane
#     local_point = plane.toLocalCoords(cq.Vector(*axis_point))

#     # Axis direction in local coordinates
#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
#     local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

#     p0 = (local_point.x, local_point.y)
#     p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

#     return profile_wp.revolve(angle, p0, p1)


# def revolve_profile(
#     profile_wp: cq.Sketch,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     wire = cq.Wire.assembleEdges(profile_wp._edges)
#     wp = cq.Workplane("XY")
#     wp.ctx.pendingWires.append(wire)

#     plane = wp.plane
#     local_point = plane.toLocalCoords(cq.Vector(*axis_point))

#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
#     local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

#     p0 = (local_point.x, local_point.y)
#     p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

#     return wp.revolve(angle, p0, p1)




# FUNCTIONA CON 12345
# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     wire = profile_wp.val()
#     if not isinstance(wire, cq.Wire):
#         raise TypeError("Expected a closed Wire on the stack (did you call .close()?)")

#     plane = profile_wp.plane
#     local_point = plane.toLocalCoords(cq.Vector(*axis_point))

#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
#     local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

#     p0 = (local_point.x, local_point.y)
#     p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

#     return profile_wp.revolve(angle, p0, p1)



# def revolve_profile(
#     profile_wp: cq.Workplane,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     sketch = profile_wp.val()
#     if not isinstance(sketch, cq.Sketch):
#         raise TypeError("Expected a cq.Sketch on the Workplane stack")

#     wire = cq.Wire.assembleEdges(sketch._edges)
#     profile_wp = cq.Workplane(profile_wp.plane).add(wire)

#     plane = profile_wp.plane
#     local_point = plane.toLocalCoords(cq.Vector(*axis_point))

#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
#     local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

#     p0 = (local_point.x, local_point.y)
#     p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

#     return profile_wp.revolve(angle, p0, p1)


# def revolve_profile(
#     profile: cq.Sketch,
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
#     plane: PlaneName = "XY",
# ) -> cq.Workplane:
#     if not isinstance(profile, cq.Sketch):
#         raise TypeError("Expected a cq.Sketch")

#     wire = cq.Wire.assembleEdges(profile._edges)
#     wp = cq.Workplane(plane)
#     wp.ctx.pendingWires.append(wire)

#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")

#     p0 = axis_point[:2]
#     p1 = (axis_point[0] + _DIRS[axis].x, axis_point[1] + _DIRS[axis].y)

#     return wp.revolve(angle, p0, p1)



# def revolve_profile(
#     profile: Tuple[cq.Sketch, PlaneName],
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     sketch, plane = profile
#     if not isinstance(sketch, cq.Sketch):
#         raise TypeError("Expected a cq.Sketch")

#     wire = cq.Wire.assembleEdges(sketch._edges)
#     wp = cq.Workplane(plane)
#     wp.ctx.pendingWires.append(wire)

#     _DIRS = {
#         "X": cq.Vector(1, 0, 0),
#         "Y": cq.Vector(0, 1, 0),
#         "Z": cq.Vector(0, 0, 1),
#     }
#     if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")

#     p0 = axis_point[:2]
#     p1 = (axis_point[0] + _DIRS[axis].x, axis_point[1] + _DIRS[axis].y)

#     return wp.revolve(angle, p0, p1)











# def revolve_profile(
#     profile: Tuple[cq.Sketch, PlaneName],
#     angle: float = 360.0,
#     axis: AxisName = "Z",
#     axis_point: Tuple[float, float, float] = (0, 0, 0),
# ) -> cq.Workplane:
#     sketch, plane = profile
#     if not isinstance(sketch, cq.Sketch):
#         raise TypeError("Expected a cq.Sketch")

#     wire = cq.Wire.assembleEdges(sketch._edges)
#     wp = cq.Workplane(plane)
#     wp.ctx.pendingWires.append(wire)

#     # Map axis to 2D direction within the given plane
#     _AXIS_2D = {
#         "XY": {"X": (1,0), "Y": (0,1), "Z": None},
#         "XZ": {"X": (1,0), "Y": None,  "Z": (0,1)},
#         "YZ": {"X": None,  "Y": (1,0), "Z": (0,1)},
#     }
#     dir2d = _AXIS_2D[plane][axis]
#     if dir2d is None:
#         raise ValueError(f"axis='{axis}' is not in plane='{plane}'")

#     p0 = axis_point[:2]
#     p1 = (p0[0] + dir2d[0], p0[1] + dir2d[1])

#     return wp.revolve(angle, p0, p1)










"""
def revolve_profile_global(
    profile_wp: cq.Workplane,
    angle: float = 360.0,
    axis: AxisName = "Z",
    axis_point: Tuple[float, float, float] = (0, 0, 0),
) -> cq.Workplane:
"""
    #Revolve a closed 2D profile around a global axis.
    #Works with Wire-based profiles (create_wire_workplane_straight_connections).
    #axis_point: a point the axis passes through, in global coordinates.
"""
    wire = profile_wp.val()
    if not isinstance(wire, cq.Wire):
        raise TypeError("Expected a closed Wire on the stack (did you call .close()?)")

    # Convert global axis_point to local workplane coordinates
    plane = profile_wp.plane
    local_point = plane.toLocalCoords(cq.Vector(*axis_point))

    # Axis direction in local coordinates
    _DIRS = {
        "X": cq.Vector(1, 0, 0),
        "Y": cq.Vector(0, 1, 0),
        "Z": cq.Vector(0, 0, 1),
    }
    if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
    local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

    p0 = (local_point.x, local_point.y)
    p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

    return profile_wp.revolve(angle, p0, p1)



def revolve_profile(
    profile: Tuple[cq.Sketch, PlaneName],
    angle: float = 360.0,
    axis: AxisName = "Z",
    axis_point: Tuple[float, float, float] = (0, 0, 0),
) -> cq.Workplane:
    sketch, plane = profile

    wire = cq.Wire.assembleEdges(sketch._edges)
    pts = [(v.X, v.Y) for v in wire.Vertices()]
    wp = cq.Workplane(plane).polyline(pts).close()

    return revolve_profile_global(wp, angle, axis, axis_point)
"""



"""
def revolve_profile(
    profile: Tuple[cq.Sketch, PlaneName],
    angle: float = 360.0,
    axis: AxisName = "Z",
    axis_point: Tuple[float, float, float] = (0, 0, 0),
) -> cq.Workplane:
    sketch, plane_name = profile

    wire = cq.Wire.assembleEdges(sketch._edges)
    pts = [(v.X, v.Y) for v in wire.Vertices()]
    wp = cq.Workplane(plane_name).polyline(pts).close()

    wire = wp.val()
    if not isinstance(wire, cq.Wire):
        raise TypeError("Expected a closed Wire on the stack (did you call .close()?)")

    plane = wp.plane
    local_point = plane.toLocalCoords(cq.Vector(*axis_point))

    _DIRS = {
        "X": cq.Vector(1, 0, 0),
        "Y": cq.Vector(0, 1, 0),
        "Z": cq.Vector(0, 0, 1),
    }
    if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
    local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

    p0 = (local_point.x, local_point.y)
    p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

    return wp.revolve(angle, p0, p1)
"""
"""
def revolve_profile(
    profile: cq.Workplane,
    angle: float = 360.0,
    axis: AxisName = "Z",
    axis_point: Tuple[float, float, float] = (0, 0, 0),
) -> cq.Workplane:
    sketch = profile.val()
    if not isinstance(sketch, cq.Sketch):
        raise TypeError("Expected a cq.Workplane with a cq.Sketch on the stack")

    # Extract face from the placed sketch via the Compound
    faces = sketch._faces.Faces()
    if not faces:
        raise ValueError("Sketch has no faces — was it closed?")
    face = faces[0]

    plane_name = getattr(profile, '_plane_name', 'XY')
    wire = face.outerWire()
    pts = [(v.X, v.Y) for v in wire.Vertices()]
    wp = cq.Workplane(plane_name).polyline(pts).close()

    plane = wp.plane
    local_point = plane.toLocalCoords(cq.Vector(*axis_point))

    _DIRS = {
        "X": cq.Vector(1, 0, 0),
        "Y": cq.Vector(0, 1, 0),
        "Z": cq.Vector(0, 0, 1),
    }
    if axis not in _DIRS: raise ValueError("axis must be 'X', 'Y', or 'Z'")
    local_dir = plane.toLocalCoords(plane.origin + _DIRS[axis]) - plane.toLocalCoords(plane.origin)

    p0 = (local_point.x, local_point.y)
    p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)

    return wp.revolve(angle, p0, p1)
"""

# ------------------------------ SWEEP ----------------------------------

"""
def sweep_profile(
    profile_wp: cq.Workplane,
    path_func: Callable[[float], Tuple[float, float, float]],
    num_path_points: int = 50,
) -> cq.Workplane:
    # Sweep a closed 2D profile along a 3D path
    
    path_start = path_func(0)
    if not np.allclose(path_start, [0, 0, 0], atol=1e-6):
        raise ValueError(f"Path must start at origin (0,0,0), got {path_start}")
    
    t_vals = np.linspace(0, 1, num_path_points)
    path_pts = [cq.Vector(*path_func(t)) for t in t_vals]
    path_edge = cq.Edge.makeSpline(path_pts)
    path_wire = cq.Wire.assembleEdges([path_edge])
    
    logger.debug(f"Path spline created with {len(path_pts)} points")
    return profile_wp.sweep(path_wire)
"""