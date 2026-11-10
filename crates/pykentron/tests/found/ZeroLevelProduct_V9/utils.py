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


# angles in degrees; X=roll, Y=pitch, Z=yaw, all through the shape's own centre
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

    # Use whichever two local coordinates carry the axis direction
    if abs(local_dir.z) > abs(local_dir.x) and abs(local_dir.z) > abs(local_dir.y):
        # axis direction is mostly in local Z — use x and z instead
        p0 = (local_point.x, local_point.z)
        p1 = (local_point.x + local_dir.x, local_point.z + local_dir.z)
    else:
        p0 = (local_point.x, local_point.y)
        p1 = (local_point.x + local_dir.x, local_point.y + local_dir.y)



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


def insert_into(base: cq.Workplane, insert: cq.Workplane) -> cq.Workplane:
    cutter = getattr(insert, "_outer", insert)  # use outer solid if available, else insert as-is
    return base.cut(cutter).union(insert)