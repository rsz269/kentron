"""
Top plate / upper head plate builder for reactor vessel models.

Creates a flat circular plate with arbitrary hole patterns (for IHX penetrations,
control rod drives, instrument nozzles, etc.).

Positioning follows the same convention as build_solid and zz_utils:
  - plate is built centred at the origin
  - rotation_angles (roll, pitch, yaw) in degrees are applied first
  - center_coords (or center_coords_pol) moves the centroid to the final position

Hole layouts:
  - symmetric:          evenly-spaced around a circle of given radius
  - custom_angles:      explicit angles at a fixed radius
  - explicit_positions: fully explicit XY positions (offsets from plate centre)

Example
-------
>>> plate = create_top_plate(
...     plate_outer_d   = 92.0,
...     plate_thickness = 3.0,
...     center_coords   = (0, 0, 91.5),   # centroid at z = straight_h + thickness/2
...     hole_groups=[
...         dict(hole_diameter=10.0, layout="explicit_positions", positions=[(0.0, 0.0)]),
...         dict(hole_diameter=6.5,  layout="symmetric",           count=6, placement_radius=31.0),
...         dict(hole_diameter=2.0,  layout="custom_angles",       angles_deg=[30.0, 150.0, 270.0], placement_radius=20.0),
...     ],
... )
"""

from __future__ import annotations

import math
from typing import Any, Tuple

import cadquery as cq

from utils import rotate_rpy_about_self_global_axes, move_center_to, convert_polar_to_cartesian


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _hole_centers(group: dict[str, Any]) -> list[tuple[float, float]]:
    """
    Resolve a hole-group spec into (x, y) offsets from the plate centre.
    The plate is always built at the origin, so these are absolute positions
    in the pre-rotation frame.

    Keys per layout:
      symmetric        : count, placement_radius  [+ start_angle_deg, default 0]
      custom_angles    : angles_deg (list[float]), placement_radius
      explicit_positions : positions (list[tuple[float, float]])
    """
    layout = group.get("layout", "symmetric")

    if layout == "symmetric":
        n     = int(group["count"])
        r     = float(group["placement_radius"])
        start = float(group.get("start_angle_deg", 0.0))
        return [
            (
                r * math.cos(math.radians(start + 360.0 * i / n)),
                r * math.sin(math.radians(start + 360.0 * i / n)),
            )
            for i in range(n)
        ]

    elif layout == "custom_angles":
        r = float(group["placement_radius"])
        return [
            (
                r * math.cos(math.radians(a)),
                r * math.sin(math.radians(a)),
            )
            for a in group["angles_deg"]
        ]

    elif layout == "explicit_positions":
        return [(float(x), float(y)) for x, y in group["positions"]]

    else:
        raise ValueError(
            f"Unknown hole layout '{layout}'. "
            "Use 'symmetric', 'custom_angles', or 'explicit_positions'."
        )


def _make_cutter(hole_d: float, plate_thickness: float, x: float, y: float) -> cq.Workplane:
    """
    Return a solid cylinder centred at (x, y, 0) — through the full plate
    thickness which is centred at z=0.  Used as a Boolean cutter.
    """
    EPS = 1e-3
    cut_h = plate_thickness + 2.0 * EPS
    # # ---- Removing the following because it creates a circular import with build_solid()
    # cutter, _ = build_solid(
    #     "primitive",
    #     {"obj_type": "cylinder", "height": cut_h, "radius": hole_d / 2.0},
    #     center_coords=(x, y, 0.0),
    # )
    # return cutter
    return cq.Workplane("XY").cylinder(cut_h, hole_d / 2.0).translate((x, y, 0.0))

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def create_top_plate(
    plate_outer_d: float,
    plate_thickness: float,
    hole_groups: list[dict[str, Any]] | None = None,
    *,
    center_coords:     Tuple[float, float, float] | None = None,
    center_coords_pol: Tuple[float, float, float] | None = None,
    rotation_angles:   Tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> cq.Workplane:
    """
    Build a flat circular plate with punched-through holes.

    The plate is built centred at the origin, then rotation and translation
    are applied — identical to build_solid positioning.

    Parameters
    ----------
    plate_outer_d : float
        Outer diameter of the plate.
    plate_thickness : float
        Thickness (height) of the plate.
    hole_groups : list[dict], optional
        Each element describes one group of identically-sized holes.  Keys:

          hole_diameter  (float) – diameter of every hole in this group
          layout         (str)   – 'symmetric' | 'custom_angles' | 'explicit_positions'

          For 'symmetric':
            count            (int)   – number of holes
            placement_radius (float) – radial distance from plate centre
            start_angle_deg  (float) – first hole angle in degrees (default 0)

          For 'custom_angles':
            angles_deg       (list[float]) – explicit angles in degrees
            placement_radius (float)

          For 'explicit_positions':
            positions (list[tuple[float, float]]) – XY offsets from plate centre

    center_coords : tuple[float, float, float], optional
        (x, y, z) of the plate centroid after rotation.
    center_coords_pol : tuple[float, float, float], optional
        (r, theta_rad, z) polar alternative to center_coords.
    rotation_angles : tuple[float, float, float]
        (roll, pitch, yaw) in degrees applied before translation.

    Returns
    -------
    cq.Workplane
        Solid plate with all holes cut through, positioned and rotated.
    """
    if plate_outer_d <= 0:
        raise ValueError("plate_outer_d must be > 0")
    if plate_thickness <= 0:
        raise ValueError("plate_thickness must be > 0")

    # ------------------------------------------------------------------ #
    # 1.  Build base plate centred at origin                              #
    # ------------------------------------------------------------------ #
    # # ---- Removing the following because it creates a circular import with build_solid()
    # plate_solid, _ = build_solid(
    #     "primitive",
    #     {
    #         "obj_type": "cylinder",
    #         "height":   plate_thickness,
    #         "radius":   plate_outer_d / 2.0,
    #     },
    # )
    plate_solid = cq.Workplane("XY").cylinder(plate_thickness, plate_outer_d / 2.0)

    # ------------------------------------------------------------------ #
    # 2.  Cut hole groups                                                  #
    # ------------------------------------------------------------------ #
    for g_idx, group in enumerate(hole_groups or []):
        hole_d = float(group["hole_diameter"])
        if hole_d <= 0:
            raise ValueError(f"hole_groups[{g_idx}]: hole_diameter must be > 0")

        for hx, hy in _hole_centers(group):
            dist = math.hypot(hx, hy)
            if dist + hole_d / 2.0 > plate_outer_d / 2.0:
                import warnings
                warnings.warn(
                    f"hole_groups[{g_idx}]: hole at ({hx:.2f}, {hy:.2f}) extends "
                    "outside the plate boundary — check placement_radius / centers_xy.",
                    stacklevel=2,
                )
            plate_solid = plate_solid.cut(
                _make_cutter(hole_d, plate_thickness, hx, hy)
            )

    # ------------------------------------------------------------------ #
    # 3.  Rotate then translate — same pattern as build_solid             #
    # ------------------------------------------------------------------ #
    if center_coords_pol is not None:
        center_coords = convert_polar_to_cartesian(*center_coords_pol)

    roll, pitch, yaw = rotation_angles
    plate_solid = rotate_rpy_about_self_global_axes(plate_solid, roll, pitch, yaw) # type: ignore

    if center_coords is not None:
        plate_solid = move_center_to(plate_solid, center_coords)

    return plate_solid


# ---------------------------------------------------------------------------
# Usage examples
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import cadquery as cq
    from ocp_vscode import show

    # ------------------------------------------------------------------
    # Example A: plate sitting on top of a vessel (inner_d=90, wall_t=1,
    # straight_h=90).  Centroid at z = 90 + 3/2 = 91.5
    # ------------------------------------------------------------------
    plate_A = create_top_plate(
        plate_outer_d   = 92.0,
        plate_thickness = 3.0,
        center_coords   = (0.0, 0.0, 91.5),
        hole_groups=[
            dict(hole_diameter=10.0, layout="explicit_positions", positions=[(0.0, 0.0)]),
            dict(hole_diameter=6.5,  layout="symmetric",          count=6, placement_radius=31.0),
            dict(hole_diameter=2.0,  layout="custom_angles",      angles_deg=[30.0, 150.0, 270.0], placement_radius=20.0),
            dict(hole_diameter=4.0,  layout="explicit_positions", positions=[(38.0, 0.0), (-38.0, 0.0)]),
        ],
    )

    # ------------------------------------------------------------------
    # Example B: 4-fold symmetric holes, tilted 15 deg in pitch
    # ------------------------------------------------------------------
    plate_B = create_top_plate(
        plate_outer_d   = 92.0,
        plate_thickness = 2.0,
        center_coords   = (0.0, 0.0, 0.0),
        rotation_angles = (0.0, 15.0, 0.0),
        hole_groups=[
            dict(hole_diameter=8.0, layout="symmetric", count=4, placement_radius=40.0, start_angle_deg=45.0),
        ],
    )

    assembly = cq.Assembly()
    assembly.add(plate_A,                          name="plate_A")
    assembly.add(plate_B.translate((110.0, 0, 0)), name="plate_B")

    show(assembly)