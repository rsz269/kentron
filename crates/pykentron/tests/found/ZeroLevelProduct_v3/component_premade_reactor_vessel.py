"""
Reactor pressure vessel (RPV) builder.

Constructs a hollow cylindrical vessel with a choice of bottom head geometry
and an optional top closure plate with arbitrary hole patterns.

The vessel shell always runs from z=0 to z=straight_h, so the top face is
exactly at z=straight_h — no bounding box lookup needed. The top plate
centroid is placed at z=straight_h + plate_thickness/2.

Returns (vessel, top_plate) as a tuple. top_plate is None if
top_plate_thickness is not provided.

Example
-------
>>> vessel, top_plate = create_reactor_vessel(
...     inner_d    = 4.72,
...     wall_t     = 0.04,
...     straight_h = 5.5,
...     bottom_head_type   = "ellipsoidal",
...     bottom_head_params = {"head_depth": 1.0},
...     top_plate_thickness = 0.1,
...     top_plate_hole_groups=[
...         dict(hole_diameter=0.52, layout="custom_angles",
...              angles_deg=[0.0, 90.0, 180.0, 270.0], placement_radius=1.7),
...     ],
... )
"""

from __future__ import annotations

import math
from typing import Any

import cadquery as cq

from component_premade_top_plate import create_top_plate


# ---------------------------------------------------------------------------
# Bottom head builders
# Each returns a solid with its rim in the XY plane at z=0,
# extending downward (z < 0).
# ---------------------------------------------------------------------------

def _head_flat(od: float, plate_t: float, **_) -> cq.Workplane:
    """Flat circular end plate, z in [-plate_t, 0]."""
    return cq.Workplane("XY").circle(od / 2).extrude(-plate_t)


def _head_hemispherical(od: float, **_) -> cq.Workplane:
    """Hemispherical bottom head, z in [-od/2, 0]."""
    r = od / 2
    sphere = cq.Workplane("XY").sphere(r)
    cutter = cq.Workplane("XY").box(10 * od, 10 * od, 10 * od).translate((0, 0, -5 * od))
    return sphere.intersect(cutter)


def _head_ellipsoidal(od: float, head_depth: float, n: int = 40, **_) -> cq.Workplane:
    """
    Ellipsoidal bottom head.
    head_depth is the dish depth (for 2:1 ellipsoidal use head_depth = od/4).
    """
    r = od / 2
    pts = [
        (r * math.cos(t), -head_depth * math.sin(t))
        for t in [i * (math.pi / 2) / n for i in range(n + 1)]
    ]
    prof = (
        cq.Workplane("XZ")
        .moveTo(*pts[0])
        .spline(pts[1:], includeCurrent=True)
        .lineTo(0, 0)
        .close()
    )
    return prof.revolve(360)


def _head_torispherical(
    od: float,
    Rc: float | None = None,
    rk: float | None = None,
    n_crown: int = 80,
    n_knuckle: int = 40,
    **_,
) -> cq.Workplane:
    """
    Torispherical (flanged & dished) bottom head.
    Rc: crown radius (defaults to od).
    rk: knuckle radius (defaults to 0.06 * od).
    """
    r = od / 2
    Rc = od if Rc is None else Rc
    rk = 0.06 * od if rk is None else rk

    if rk <= 0:          raise ValueError("rk must be > 0")
    if Rc <= r:          raise ValueError("Rc must be > od/2")
    if Rc <= rk:         raise ValueError("Rc must be > rk")

    xk, zk = r - rk, 0.0
    d = Rc - rk
    rad = d * d - xk * xk
    if rad <= 0:
        raise ValueError("Infeasible Rc/rk for this od (try larger Rc or smaller rk)")
    zc = math.sqrt(rad)

    dx, dz = xk - 0.0, zk - zc
    L = math.hypot(dx, dz)
    ux, uz = dx / L, dz / L
    xt = Rc * ux
    zt = zc + Rc * uz

    a0 = math.atan2(0.0 - zk, r - xk)
    a1 = math.atan2(zt - zk, xt - xk)
    if a1 > a0:
        a1 -= 2 * math.pi
    knuckle_pts = [
        (xk + rk * math.cos(a), zk + rk * math.sin(a))
        for a in [a0 + i * (a1 - a0) / n_knuckle for i in range(n_knuckle + 1)]
    ]

    b0 = math.atan2(zt - zc, xt)
    b1 = -math.pi / 2
    crown_pts = [
        (Rc * math.cos(b), zc + Rc * math.sin(b))
        for b in [b0 + i * (b1 - b0) / n_crown for i in range(1, n_crown + 1)]
    ]

    pts = knuckle_pts + crown_pts
    prof = (
        cq.Workplane("XZ")
        .moveTo(*pts[0])
        .spline(pts[1:], includeCurrent=True)
        .lineTo(0, 0)
        .close()
    )
    return prof.revolve(360)


_HEAD_BUILDERS = {
    "flat":          _head_flat,
    "hemispherical": _head_hemispherical,
    "ellipsoidal":   _head_ellipsoidal,
    "torispherical": _head_torispherical,
}


def _build_outer_head(od: float, head_type: str, params: dict) -> cq.Workplane:
    """Build a bottom head solid at z=0 extending downward."""
    if head_type not in _HEAD_BUILDERS:
        raise ValueError(
            f"Unknown head type '{head_type}'. "
            f"Choose from: {list(_HEAD_BUILDERS)}"
        )
    return _HEAD_BUILDERS[head_type](od, **params)


def _build_top_head(od: float, head_type: str, params: dict, z0: float) -> cq.Workplane:
    """
    Build a top head whose rim lies at z=z0 and extends upward.
    Reuses the bottom head builders by mirroring about XY.
    """
    h = _build_outer_head(od, head_type, params)
    return h.mirror("XY").translate((0, 0, z0))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def create_reactor_vessel(
    inner_d: float,
    wall_t: float,
    straight_h: float,
    *,
    bottom_head_type:   str | None = None,
    bottom_head_params: dict | None = None,
    top_head_type:      str | None = None,
    top_head_params:    dict | None = None,
    top_plate_thickness:   float | None = None,
    top_plate_hole_groups: list[dict[str, Any]] | None = None,
) -> tuple[cq.Workplane, cq.Workplane | None]:
    """
    Build a reactor pressure vessel with optional head geometries and top plate.

    The cylindrical shell always runs from z=0 to z=straight_h. The top face
    is therefore always exactly at z=straight_h, making top plate positioning
    exact without any bounding box calculation.

    Parameters
    ----------
    inner_d : float
        Inner diameter of the cylindrical shell.
    wall_t : float
        Wall thickness. Outer diameter = inner_d + 2 * wall_t.
    straight_h : float
        Height of the straight cylindrical section.

    bottom_head_type : str, optional
        'flat' | 'hemispherical' | 'ellipsoidal' | 'torispherical'
        If None, the vessel has an open bottom.
    bottom_head_params : dict, optional
        Parameters forwarded to the bottom head builder:
          flat:          plate_t (defaults to wall_t)
          ellipsoidal:   head_depth (required)
          torispherical: Rc, rk (both optional, have defaults)
          hemispherical: no extra params needed

    top_head_type : str, optional
        Same options as bottom_head_type. If None, the vessel has an open top.
        Note: if top_plate_thickness is also given, the plate sits on top of
        the top head.
    top_head_params : dict, optional
        Same structure as bottom_head_params.

    top_plate_thickness : float, optional
        If given, a flat closure plate of this thickness is created and
        returned as the second element of the tuple. Its bottom face sits
        flush at z=straight_h (or on top of the top head if one is present).
    top_plate_hole_groups : list[dict], optional
        Passed directly to create_top_plate. See create_top_plate for the
        full hole_groups specification.

    Returns
    -------
    (vessel, top_plate) : tuple[cq.Workplane, cq.Workplane | None]
        vessel    — the RPV shell with head(s) attached
        top_plate — the closure plate, or None if top_plate_thickness is None
    """
    if inner_d <= 0:   raise ValueError("inner_d must be > 0")
    if wall_t  <= 0:   raise ValueError("wall_t must be > 0")
    if straight_h <= 0: raise ValueError("straight_h must be > 0")

    bottom_head_params = dict(bottom_head_params or {})
    top_head_params    = dict(top_head_params    or {})

    # defaults
    if bottom_head_type == "flat":
        bottom_head_params.setdefault("plate_t", wall_t)
    if top_head_type == "flat":
        top_head_params.setdefault("plate_t", wall_t)

    # validation
    if bottom_head_type == "ellipsoidal" and "head_depth" not in bottom_head_params:
        raise ValueError("bottom_head_type='ellipsoidal' requires head_depth in bottom_head_params")
    if top_head_type == "ellipsoidal" and "head_depth" not in top_head_params:
        raise ValueError("top_head_type='ellipsoidal' requires head_depth in top_head_params")

    od = inner_d + 2 * wall_t

    # ------------------------------------------------------------------ #
    # 1.  Outer shell + heads                                             #
    # ------------------------------------------------------------------ #
    outer = cq.Workplane("XY").circle(od / 2).extrude(straight_h)

    if bottom_head_type == "flat":
        t = float(bottom_head_params["plate_t"])
        outer = outer.union(cq.Workplane("XY").circle(od / 2).extrude(-t))
    elif bottom_head_type is not None:
        outer = outer.union(_build_outer_head(od, bottom_head_type, bottom_head_params))

    if top_head_type == "flat":
        t = float(top_head_params["plate_t"])
        outer = outer.union(
            cq.Workplane("XY").workplane(offset=straight_h).circle(od / 2).extrude(t)
        )
    elif top_head_type is not None:
        outer = outer.union(_build_top_head(od, top_head_type, top_head_params, straight_h))

    # ------------------------------------------------------------------ #
    # 2.  Inner bore (cutter)                                             #
    # ------------------------------------------------------------------ #
    inner = cq.Workplane("XY").circle(inner_d / 2).extrude(straight_h)

    if bottom_head_type == "ellipsoidal":
        hd = float(bottom_head_params["head_depth"])
        p = dict(bottom_head_params)
        p["head_depth"] = max(hd - wall_t, 1e-3)
        inner = inner.union(_build_outer_head(inner_d, "ellipsoidal", p))
    elif bottom_head_type not in (None, "flat"):
        inner = inner.union(_build_outer_head(inner_d, bottom_head_type, bottom_head_params))

    if top_head_type == "ellipsoidal":
        hd = float(top_head_params["head_depth"])
        p = dict(top_head_params)
        p["head_depth"] = max(hd - wall_t, 1e-3)
        inner = inner.union(_build_top_head(inner_d, "ellipsoidal", p, straight_h))
    elif top_head_type not in (None, "flat"):
        inner = inner.union(_build_top_head(inner_d, top_head_type, top_head_params, straight_h))

    # The following cut produced different objects, leading to problems in downstream operations when exporting to DAGMC
    # vessel = outer.cut(inner).clean()

    # After the boolean operations (union + cut), CadQuery/OCCT does not always
    # produce a single topological solid. Instead, it can leave the result as a
    # "compound" — multiple sub-shapes (e.g. cylinder + hemisphere) that are
    # visually merged but internally still separate entities.
    #
    # This causes problems downstream when cad_to_dagmc processes the STEP file:
    # it counts each sub-shape as a separate volume, so what looks like 1 component
    # becomes 2 volumes, causing a material tag mismatch.
    #
    # The fix below forces OCCT to perform a true topological fusion (Fuse),
    # which merges all sub-shapes into a single connected solid with no internal
    # boundaries. This ensures the STEP export contains exactly 1 volume per component.

    # The following doesn't work because The result is a Compound, not a Solid, so .Fuse() isn't available. 
    # Use CadQuery's own fusion instead. Replace the fix in reactor_vessel.py with hte other code below:

    # vessel = outer.cut(inner).clean()
    # vessel = cq.Workplane().add(
    # cq.Shape.cast(vessel.val().wrapped.Fuse(vessel.val().wrapped)) #type: ignore
    # )


    vessel = outer.cut(inner).clean()

    # Force true topological fusion into a single solid.
    # .combine() merges all sub-shapes in the compound into one connected solid.
    solids = vessel.solids().vals()
    fused = solids[0]
    for s in solids[1:]:
        fused = fused.fuse(s)  # type: ignore
    vessel = cq.Workplane().add(fused)









    # ------------------------------------------------------------------ #
    # 3.  Top plate                                                        #
    # ------------------------------------------------------------------ #
    top_plate = None
    if top_plate_thickness is not None:
        # top face is always exactly at z=straight_h — no bounding box needed
        top_plate = create_top_plate(
            plate_outer_d   = od,
            plate_thickness = top_plate_thickness,
            center_coords   = (0.0, 0.0, straight_h + top_plate_thickness / 2.0),
            hole_groups     = top_plate_hole_groups,
        )

    return vessel, top_plate


# ---------------------------------------------------------------------------
# Usage examples
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import cadquery as cq
    from ocp_vscode import show

    # ------------------------------------------------------------------
    # Example A: ellipsoidal bottom, no top head, flat top plate with
    # 4 IHX penetrations matching REACTOR_HX layout (r=1.7, 4-fold)
    # ------------------------------------------------------------------
    vessel_A, plate_A = create_reactor_vessel(
        inner_d    = 4.72,
        wall_t     = 0.04,
        straight_h = 5.5,
        bottom_head_type   = "ellipsoidal",
        bottom_head_params = {"head_depth": 1.0},
        top_plate_thickness = 0.1,
        top_plate_hole_groups=[
            dict(
                hole_diameter    = 0.52,
                layout           = "custom_angles",
                angles_deg       = [0.0, 90.0, 180.0, 270.0],
                placement_radius = 1.7,
            ),
        ],
    )

    assembly_A = cq.Assembly()
    assembly_A.add(vessel_A, name="rpv")
    if plate_A is not None:
        assembly_A.add(plate_A, name="top_plate")
    show(assembly_A)
    import time; time.sleep(7)
    # ------------------------------------------------------------------
    # Example B: torispherical bottom, flat top head, no top plate
    # ------------------------------------------------------------------
    vessel_B, _ = create_reactor_vessel(
        inner_d    = 4.72,
        wall_t     = 0.04,
        straight_h = 5.5,
        bottom_head_type   = "torispherical",
        bottom_head_params = {"Rc": 4.72, "rk": 0.06 * 4.72},
        top_head_type      = "flat",
    )

    show(vessel_B)