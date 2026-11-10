"""
Parametric SFR above-core structure (ACS) — ZLP V3 (simplified, solid bottom).

Three bodies unioned:

  Lower shell   [0  → z4]   bottom ring + cone + collar + neck
                             revolved around the LOCAL ORIGIN (cone axis).
                             SOLID — no inner bore.

  Closing plate [z4-cp_h → z4]
                             solid disk (outer r = top_cyl_outer_r),
                             translated by (top_cyl_offset_x, top_cyl_offset_y).

  Upper cylinder [z4 → z5]  plain solid cylinder,
                             translated by (top_cyl_offset_x, top_cyl_offset_y).

The cone-axis stays on the local origin so the lower shell (and the hex
through-hole pattern beneath it) can be aligned with the reactor core by
the assembly's center_coords. The top cylinder can be displaced sideways
via top_cyl_offset_x/y to clear surrounding components.

Stacking order (bottom → top):
  z=0     bottom face
  z1      bottom ring top          ( = bottom_ring_height )
  z2      cone top = collar bottom ( = z1 + cone_height   )
  z3      collar top = neck bottom ( = z2 + collar_height )
  z4      neck top = cylinder bottom ( = z3 + neck_height )
  z5      cylinder top

Single public function:  create_above_core_structure()
"""

from __future__ import annotations
import math
import cadquery as cq

from profile_from_straight_connections import create_profile_from_straight_connections
from utils import revolve_profile


def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    EPS = 1e-9
    out = [pts[0]]
    for p in pts[1:]:
        if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
            out.append(p)
    return out


def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
    profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
    return revolve_profile(profile, angle=360, axis="Z")


def create_above_core_structure(
    # ── Top cylinder — upper axis (0, 0) ─────────────────────────────────
    top_cyl_outer_r: float,  # outer radius of the top cylinder [m]
    top_cyl_height:  float,  # axial height of the top cylinder [m]

    # ── Neck — lower / cone axis ──────────────────────────────────────────
    neck_outer_r: float,
    neck_height:  float,

    # ── Collar — lower / cone axis ────────────────────────────────────────
    collar_outer_r: float,
    collar_height:  float,

    # ── Wall thickness (uniform) ──────────────────────────────────────────
    wall_t: float,

    # ── Cone + bottom ring — lower / cone axis ────────────────────────────
    cone_bottom_outer_r: float,
    cone_height:         float,
    bottom_ring_height:  float,

    # ── Closing plate ─────────────────────────────────────────────────────
    closing_plate_height: float,   # thickness of the plate at z4 [m]

    # ── Top cylinder offset relative to the lower-shell (cone) axis ──────
    # The lower shell (cone + collar + neck) sits on the local origin. The
    # top cylinder + closing plate are translated by this offset. Set to
    # (0, 0) for a coaxial component; nonzero to displace the top cylinder
    # sideways (e.g. to clear pumps / IHX nozzles in the assembly).
    top_cyl_offset_x: float,
    top_cyl_offset_y: float,

    # ── Optional flow holes on the cone ───────────────────────────────────
    flow_hole_groups: list[dict] | None = None,

    # ── Optional through-holes in the bottom (hex pattern: 1 + 6) ─────────
    # dict with keys:
    #   "through_d":    through-hole diameter [m]                   (e.g. 0.080)
    #   "counter_d":    counterbore diameter   [m]                  (e.g. 0.142)
    #   "counter_depth": counterbore depth from top of plate [m]    (default = closing_plate_height)
    #   "pitch":        center-to-center spacing of the hex ring [m] (default 0.300)
    bottom_holes: dict | None = None,

    # ── Global position ───────────────────────────────────────────────────
    z_bottom: float = 0.0,

) -> cq.Workplane:
    """
    Build the above-core structure as the union of three solids:

      1. Lower shell [0 → z4] — bottom ring + cone + collar + neck.
         Revolved around the local origin (cone axis). Fully solid.

      2. Closing plate [z4 - closing_plate_height → z4]
         Solid disk (outer r = top_cyl_outer_r), translated by
         (top_cyl_offset_x, top_cyl_offset_y).

      3. Upper cylinder [z4 → z5]
         Plain solid cylinder, translated by
         (top_cyl_offset_x, top_cyl_offset_y).
    """
    # ── Validate ─────────────────────────────────────────────────────────
    checks = [
        (collar_outer_r >= neck_outer_r,     "collar_outer_r >= neck_outer_r"),
        (neck_outer_r > wall_t,              "neck_outer_r > wall_t"),
        (cone_bottom_outer_r > neck_outer_r, "cone_bottom_outer_r > neck_outer_r"),
        (closing_plate_height > 0,           "closing_plate_height > 0"),
    ]
    for ok, msg in checks:
        if not ok:
            raise ValueError(f"Validation failed: {msg}")
    for name, val in [
        ("top_cyl_height", top_cyl_height),
        ("top_cyl_outer_r", top_cyl_outer_r),
        ("neck_height", neck_height), ("collar_height", collar_height),
        ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
        ("wall_t", wall_t),
    ]:
        if val <= 0:
            raise ValueError(f"{name} must be > 0")

    # ── Z levels ──────────────────────────────────────────────────────────
    z1 = bottom_ring_height
    z2 = z1 + cone_height
    z3 = z2 + collar_height
    z4 = z3 + neck_height        # ← axis split
    z5 = z4 + top_cyl_height

    # ── 1. Lower shell (solid — no inner bore) ────────────────────────────
    # The lower shell sits on the LOCAL ORIGIN. The top cylinder + closing
    # plate are the ones that carry top_cyl_offset_x/y. This way the cone /
    # collar / cone-axis hole pattern stay on the assembly's central axis
    # (and so can be aligned with the reactor core), while the top cylinder
    # gets displaced sideways.
    lower_pts = [
        # outer (upward)
        (cone_bottom_outer_r, 0),
        (cone_bottom_outer_r, z1),
        (neck_outer_r,        z2),
        (collar_outer_r,      z2),
        (collar_outer_r,      z3),
        (neck_outer_r,        z3),
        (neck_outer_r,        z4),
        # close along the axis back down to z=0
        (0,                   z4),
        (0,                   0),
    ]
    lower_solid = _revolve_closed(lower_pts)
    # (no translation — lower shell stays on the origin)

    # ── 2. Closing plate ──────────────────────────────────────────────────
    # Full disk; translated to the top-cylinder offset.
    closing_plate = (
        cq.Workplane("XY")
        .workplane(offset=z4 - closing_plate_height)
        .circle(top_cyl_outer_r)
        .extrude(closing_plate_height)
    )
    if top_cyl_offset_x != 0.0 or top_cyl_offset_y != 0.0:
        closing_plate = closing_plate.translate((top_cyl_offset_x, top_cyl_offset_y, 0))

    # ── 3. Upper cylinder ─────────────────────────────────────────────────
    upper_solid = (
        cq.Workplane("XY")
        .workplane(offset=z4)
        .circle(top_cyl_outer_r)
        .extrude(top_cyl_height)
    )
    if top_cyl_offset_x != 0.0 or top_cyl_offset_y != 0.0:
        upper_solid = upper_solid.translate((top_cyl_offset_x, top_cyl_offset_y, 0))

    # ── 4. Union ──────────────────────────────────────────────────────────
    solid = upper_solid.union(closing_plate).union(lower_solid)

    # ── 5. Flow holes (positions on the lower-shell / cone axis = origin) ─
    if flow_hole_groups:
        for group in flow_hole_groups:
            hole_r  = float(group["hole_r"])
            z_c     = float(group["z_center"])
            n       = int(group["n_holes"])
            start_a = float(group.get("start_angle_deg", 0.0))
            if z_c < 0 or z_c > z5:
                raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
            if z_c <= z1:
                r_mid = cone_bottom_outer_r - wall_t / 2.0
            elif z_c <= z2:
                frac  = (z_c - z1) / (z2 - z1)
                r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
                r_mid = r_out - wall_t / 2.0
            else:
                r_mid = neck_outer_r - wall_t / 2.0
            for i in range(n):
                a  = math.radians(start_a + 360.0 * i / n)
                hx = r_mid * math.cos(a)
                hy = r_mid * math.sin(a)
                cutter = (
                    cq.Workplane("XY").workplane(offset=z_c)
                    .circle(hole_r).extrude(wall_t * 2)
                    .translate((hx, hy, -wall_t))
                )
                solid = solid.cut(cutter)

    # ── 5b. Bottom hex through-holes + counterbores ──────────────────────
    # Pattern is centered on the cone axis (= local origin), so the holes
    # are radially aligned with the lower shell — and through the reactor
    # core when the assembly places the ACS on the reactor centerline.
    if bottom_holes:
        through_d     = float(bottom_holes["through_d"])
        counter_d     = float(bottom_holes["counter_d"])
        counter_depth = float(bottom_holes.get("counter_depth", closing_plate_height))
        pitch         = float(bottom_holes.get("pitch", 0.300))

        if counter_d <= through_d:
            raise ValueError("counter_d must be greater than through_d")
        if counter_depth <= 0 or counter_depth >= z5:
            raise ValueError(f"counter_depth out of range: {counter_depth}")

        through_r = through_d / 2.0
        counter_r = counter_d / 2.0

        # Hex pattern: center + 6 around it
        centers = [(0.0, 0.0)]
        for i in range(6):
            a = math.radians(60.0 * i)
            centers.append((pitch * math.cos(a), pitch * math.sin(a)))

        EPS = 1e-4
        for hx, hy in centers:
            # Through-hole: from below the bottom up through the closing plate
            # (and through the top cylinder if it happens to overlap radially).
            through_cutter = (
                cq.Workplane("XY")
                .workplane(offset=-EPS)
                .circle(through_r)
                .extrude(z5 + 2 * EPS)
                .translate((hx, hy, 0))
            )
            solid = solid.cut(through_cutter)

            # Counterbore: pocket from top of closing plate (z4) downward.
            counter_cutter = (
                cq.Workplane("XY")
                .workplane(offset=z4 - counter_depth)
                .circle(counter_r)
                .extrude(counter_depth + EPS)
                .translate((hx, hy, 0))
            )
            solid = solid.cut(counter_cutter)

    # ── 6. Final z translation ────────────────────────────────────────────
    if z_bottom != 0.0:
        solid = solid.translate((0, 0, z_bottom))

    return solid


if __name__ == "__main__":
    from ocp_vscode import show

    acs = create_above_core_structure(
        top_cyl_outer_r      = 1.843,
        top_cyl_height       = 1.008,
        neck_outer_r         = 1.1085,
        neck_height          = 0.569,
        collar_outer_r       = 1.1085,
        collar_height        = 0.092,
        wall_t               = 0.025,
        cone_height          = 2.429,
        cone_bottom_outer_r  = 1.403,
        bottom_ring_height   = 0.498,
        closing_plate_height = 0.050,
        top_cyl_offset_x     = 0.6056,
        top_cyl_offset_y     = 0.0,
        bottom_holes = {
            "through_d":     0.080,   # Ø80 mm through-hole
            "counter_d":     0.142,   # Ø142 mm counterbore
            "counter_depth": 0.050,   # depth of counterbore from top of plate
            "pitch":         0.300,   # center-to-center of hex ring
        },
    )
    show(acs)
























































































# # 19052026 15:20

# """
# Parametric SFR above-core structure (ACS) — ZLP V3 (simplified, solid bottom).

# Three bodies unioned:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              SOLID — no inner bore.

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_cyl_outer_r, origin).

#   Upper cylinder [z4 → z5]  plain solid cylinder around (0, 0).

# Stacking order (bottom → top):
#   z=0     bottom face
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = cylinder bottom ( = z3 + neck_height ) ← axis split
#   z5      cylinder top

# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Top cylinder — upper axis (0, 0) ─────────────────────────────────
#     top_cyl_outer_r: float,  # outer radius of the top cylinder [m]
#     top_cyl_height:  float,  # axial height of the top cylinder [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Optional through-holes in the bottom (hex pattern: 1 + 6) ─────────
#     # dict with keys:
#     #   "through_d":    through-hole diameter [m]                   (e.g. 0.080)
#     #   "counter_d":    counterbore diameter   [m]                  (e.g. 0.142)
#     #   "counter_depth": counterbore depth from top of plate [m]    (default = closing_plate_height)
#     #   "pitch":        center-to-center spacing of the hex ring [m] (default 0.300)
#     bottom_holes: dict | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Fully solid (no inner bore).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_cyl_outer_r, centered at origin).

#       3. Upper cylinder [z4 → z5]
#          Plain solid cylinder around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (collar_outer_r >= neck_outer_r,     "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,              "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r, "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,           "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_cyl_height", top_cyl_height),
#         ("top_cyl_outer_r", top_cyl_outer_r),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1 = bottom_ring_height
#     z2 = z1 + cone_height
#     z3 = z2 + collar_height
#     z4 = z3 + neck_height        # ← axis split
#     z5 = z4 + top_cyl_height

#     # ── 1. Lower shell (solid — no inner bore) ────────────────────────────
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # close along the axis back down to z=0
#         (0,                   z4),
#         (0,                   0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin. Bottom is solid, so no bore cut is needed.
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_cyl_outer_r)
#         .extrude(closing_plate_height)
#     )

#     # ── 3. Upper cylinder ─────────────────────────────────────────────────
#     upper_solid = (
#         cq.Workplane("XY")
#         .workplane(offset=z4)
#         .circle(top_cyl_outer_r)
#         .extrude(top_cyl_height)
#     )

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 5b. Bottom hex through-holes + counterbores ──────────────────────
#     if bottom_holes:
#         through_d     = float(bottom_holes["through_d"])
#         counter_d     = float(bottom_holes["counter_d"])
#         counter_depth = float(bottom_holes.get("counter_depth", closing_plate_height))
#         pitch         = float(bottom_holes.get("pitch", 0.300))

#         if counter_d <= through_d:
#             raise ValueError("counter_d must be greater than through_d")
#         if counter_depth <= 0 or counter_depth >= z5:
#             raise ValueError(f"counter_depth out of range: {counter_depth}")

#         through_r = through_d / 2.0
#         counter_r = counter_d / 2.0

#         # Hex pattern: center + 6 around it
#         centers = [(0.0, 0.0)]
#         for i in range(6):
#             a = math.radians(60.0 * i)
#             centers.append((pitch * math.cos(a), pitch * math.sin(a)))

#         # Translate hex centers into world coordinates (cone-axis frame)
#         EPS = 1e-4
#         for dx, dy in centers:
#             hx = cx + dx
#             hy = cy + dy

#             # Through-hole: from below the bottom up through the top cylinder
#             through_cutter = (
#                 cq.Workplane("XY")
#                 .workplane(offset=-EPS)
#                 .circle(through_r)
#                 .extrude(z5 + 2 * EPS)
#                 .translate((hx, hy, 0))
#             )
#             solid = solid.cut(through_cutter)

#             # Counterbore: pocket from top of plate (z4) downward by counter_depth
#             counter_cutter = (
#                 cq.Workplane("XY")
#                 .workplane(offset=z4 - counter_depth)
#                 .circle(counter_r)
#                 .extrude(counter_depth + EPS)
#                 .translate((hx, hy, 0))
#             )
#             solid = solid.cut(counter_cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_cyl_outer_r      = 1.843,
#         top_cyl_height       = 1.008,
#         neck_outer_r         = 1.1085,
#         neck_height          = 0.569,
#         collar_outer_r       = 1.1085,
#         collar_height        = 0.092,
#         wall_t               = 0.025,
#         cone_height          = 2.429,
#         cone_bottom_outer_r  = 1.403,
#         bottom_ring_height   = 0.498,
#         closing_plate_height = 0.050,
#         cone_axis_offset_x   = 0.6056,
#         cone_axis_offset_y   = 0.0,
#         bottom_holes = {
#             "through_d":     0.080,   # Ø80 mm through-hole
#             "counter_d":     0.142,   # Ø142 mm counterbore
#             "counter_depth": 0.050,   # depth of counterbore from top of plate
#             "pitch":         0.300,   # center-to-center of hex ring
#         },
#     )
#     show(acs)








































































































# """
# Parametric SFR above-core structure (ACS) — ZLP V3 (simplified, solid bottom).

# Three bodies unioned:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              SOLID — no inner bore.

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_cyl_outer_r, origin).

#   Upper cylinder [z4 → z5]  plain solid cylinder around (0, 0).

# Stacking order (bottom → top):
#   z=0     bottom face
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = cylinder bottom ( = z3 + neck_height ) ← axis split
#   z5      cylinder top

# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Top cylinder — upper axis (0, 0) ─────────────────────────────────
#     top_cyl_outer_r: float,  # outer radius of the top cylinder [m]
#     top_cyl_height:  float,  # axial height of the top cylinder [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Optional through-holes in the bottom (hex pattern: 1 + 6) ─────────
#     # dict with keys:
#     #   "through_d":    through-hole diameter [m]                   (e.g. 0.080)
#     #   "counter_d":    counterbore diameter   [m]                  (e.g. 0.142)
#     #   "counter_depth": counterbore depth from top of plate [m]    (default = closing_plate_height)
#     #   "pitch":        center-to-center spacing of the hex ring [m] (default 0.300)
#     bottom_holes: dict | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Fully solid (no inner bore).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_cyl_outer_r, centered at origin).

#       3. Upper cylinder [z4 → z5]
#          Plain solid cylinder around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (collar_outer_r >= neck_outer_r,     "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,              "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r, "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,           "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_cyl_height", top_cyl_height),
#         ("top_cyl_outer_r", top_cyl_outer_r),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1 = bottom_ring_height
#     z2 = z1 + cone_height
#     z3 = z2 + collar_height
#     z4 = z3 + neck_height        # ← axis split
#     z5 = z4 + top_cyl_height

#     # ── 1. Lower shell (solid — no inner bore) ────────────────────────────
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # close along the axis back down to z=0
#         (0,                   z4),
#         (0,                   0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin. Bottom is solid, so no bore cut is needed.
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_cyl_outer_r)
#         .extrude(closing_plate_height)
#     )

#     # ── 3. Upper cylinder ─────────────────────────────────────────────────
#     upper_solid = (
#         cq.Workplane("XY")
#         .workplane(offset=z4)
#         .circle(top_cyl_outer_r)
#         .extrude(top_cyl_height)
#     )

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 5b. Bottom hex through-holes + counterbores ──────────────────────
#     if bottom_holes:
#         through_d     = float(bottom_holes["through_d"])
#         counter_d     = float(bottom_holes["counter_d"])
#         counter_depth = float(bottom_holes.get("counter_depth", closing_plate_height))
#         pitch         = float(bottom_holes.get("pitch", 0.300))

#         if counter_d <= through_d:
#             raise ValueError("counter_d must be greater than through_d")
#         if counter_depth <= 0 or counter_depth >= z5:
#             raise ValueError(f"counter_depth out of range: {counter_depth}")

#         through_r = through_d / 2.0
#         counter_r = counter_d / 2.0

#         # Hex pattern: center + 6 around it
#         centers = [(0.0, 0.0)]
#         for i in range(6):
#             a = math.radians(60.0 * i)
#             centers.append((pitch * math.cos(a), pitch * math.sin(a)))

#         # Translate hex centers into world coordinates (cone-axis frame)
#         EPS = 1e-4
#         for dx, dy in centers:
#             hx = cx + dx
#             hy = cy + dy

#             # Through-hole: from below the bottom up to z4 (top of closing plate)
#             through_cutter = (
#                 cq.Workplane("XY")
#                 .workplane(offset=-EPS)
#                 .circle(through_r)
#                 .extrude(z4 + EPS)
#                 .translate((hx, hy, 0))
#             )
#             solid = solid.cut(through_cutter)

#             # Counterbore: pocket from top of plate (z4) downward by counter_depth
#             counter_cutter = (
#                 cq.Workplane("XY")
#                 .workplane(offset=z4 - counter_depth)
#                 .circle(counter_r)
#                 .extrude(counter_depth + EPS)
#                 .translate((hx, hy, 0))
#             )
#             solid = solid.cut(counter_cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_cyl_outer_r      = 1.843,
#         top_cyl_height       = 1.008,
#         neck_outer_r         = 1.1085,
#         neck_height          = 0.569,
#         collar_outer_r       = 1.1085,
#         collar_height        = 0.092,
#         wall_t               = 0.025,
#         cone_height          = 2.429,
#         cone_bottom_outer_r  = 1.403,
#         bottom_ring_height   = 0.498,
#         closing_plate_height = 0.050,
#         cone_axis_offset_x   = 0.6056,
#         cone_axis_offset_y   = 0.0,
#         bottom_holes = {
#             "through_d":     0.080,   # Ø80 mm through-hole
#             "counter_d":     0.142,   # Ø142 mm counterbore
#             "counter_depth": 0.050,   # depth of counterbore from top of plate
#             "pitch":         0.300,   # center-to-center of hex ring
#         },
#     )
#     show(acs)






















































# """
# Parametric SFR above-core structure (ACS) — ZLP V3 (simplified, solid bottom).

# Three bodies unioned:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              SOLID — no inner bore.

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_cyl_outer_r, origin).

#   Upper cylinder [z4 → z5]  plain solid cylinder around (0, 0).

# Stacking order (bottom → top):
#   z=0     bottom face
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = cylinder bottom ( = z3 + neck_height ) ← axis split
#   z5      cylinder top

# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Top cylinder — upper axis (0, 0) ─────────────────────────────────
#     top_cyl_outer_r: float,  # outer radius of the top cylinder [m]
#     top_cyl_height:  float,  # axial height of the top cylinder [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Fully solid (no inner bore).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_cyl_outer_r, centered at origin).

#       3. Upper cylinder [z4 → z5]
#          Plain solid cylinder around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (collar_outer_r >= neck_outer_r,     "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,              "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r, "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,           "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_cyl_height", top_cyl_height),
#         ("top_cyl_outer_r", top_cyl_outer_r),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1 = bottom_ring_height
#     z2 = z1 + cone_height
#     z3 = z2 + collar_height
#     z4 = z3 + neck_height        # ← axis split
#     z5 = z4 + top_cyl_height

#     # ── 1. Lower shell (solid — no inner bore) ────────────────────────────
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # close along the axis back down to z=0
#         (0,                   z4),
#         (0,                   0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin. Bottom is solid, so no bore cut is needed.
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_cyl_outer_r)
#         .extrude(closing_plate_height)
#     )

#     # ── 3. Upper cylinder ─────────────────────────────────────────────────
#     upper_solid = (
#         cq.Workplane("XY")
#         .workplane(offset=z4)
#         .circle(top_cyl_outer_r)
#         .extrude(top_cyl_height)
#     )

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_cyl_outer_r      = 1.843,
#         top_cyl_height       = 1.008,
#         neck_outer_r         = 1.1085,
#         neck_height          = 0.569,
#         collar_outer_r       = 1.1085,
#         collar_height        = 0.092,
#         wall_t               = 0.025,
#         cone_height          = 2.429,
#         cone_bottom_outer_r  = 1.403,
#         bottom_ring_height   = 0.498,
#         closing_plate_height = 0.050,
#         cone_axis_offset_x   = 0.6056,
#         cone_axis_offset_y   = 0.0,
#     )
#     show(acs)




























































# """
# Parametric SFR above-core structure (ACS) — ZLP V3 (simplified top).

# Three bodies unioned, one bore cut:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              bottom face CLOSED (profile goes to r=0 at z=0)

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_cyl_outer_r, origin)
#                              with neck inner bore cut through at cone-axis offset

#   Upper cylinder [z4 → z5]  plain solid cylinder around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face (closed)
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = cylinder bottom ( = z3 + neck_height ) ← axis split
#   z5      cylinder top

# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Top cylinder — upper axis (0, 0) ─────────────────────────────────
#     top_cyl_outer_r: float,  # outer radius of the top cylinder [m]
#     top_cyl_height:  float,  # axial height of the top cylinder [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids minus one bore:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Bottom face is closed (inner profile sweeps to r=0 at z=0).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_cyl_outer_r, centered at origin).
#          The neck inner bore (neck_inner_r at cone-axis offset) is cut through.

#       3. Upper cylinder [z4 → z5]
#          Plain solid cylinder around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (collar_outer_r >= neck_outer_r,     "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,              "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r, "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,           "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_cyl_height", top_cyl_height),
#         ("top_cyl_outer_r", top_cyl_outer_r),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1 = bottom_ring_height
#     z2 = z1 + cone_height
#     z3 = z2 + collar_height
#     z4 = z3 + neck_height        # ← axis split
#     z5 = z4 + top_cyl_height

#     neck_inner_r = neck_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     # Inner profile descends to (0, 0) to close the bottom face.
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # inner (downward) — bore extends all the way to z=0
#         (neck_inner_r,        z4),
#         (neck_inner_r,        0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin; neck inner bore cut through at cone-axis offset.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_cyl_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper cylinder ─────────────────────────────────────────────────
#     upper_solid = (
#         cq.Workplane("XY")
#         .workplane(offset=z4)
#         .circle(top_cyl_outer_r)
#         .extrude(top_cyl_height)
#     )

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_cyl_outer_r      = 1.843,
#         top_cyl_height       = 1.008,
#         neck_outer_r         = 1.1085,
#         neck_height          = 0.569,
#         collar_outer_r       = 1.1085,
#         collar_height        = 0.092,
#         wall_t               = 0.025,
#         cone_height          = 2.429,
#         cone_bottom_outer_r  = 1.403,
#         bottom_ring_height   = 0.498,
#         closing_plate_height = 0.050,
#         cone_axis_offset_x   = 0.6056,
#         cone_axis_offset_y   = 0.0,
#     )
#     show(acs)































# """
# Parametric SFR above-core structure (ACS) — ZLP V3.

# Three bodies unioned, one bore cut:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              bottom face CLOSED (profile goes to r=0 at z=0)

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_flange_outer_r, origin)
#                              with neck inner bore cut through at cone-axis offset

#   Upper shell   [z4 → z5]   flange (stepped outer profile, inner bore open)
#                              revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face (closed)
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = flange bottom ( = z3 + neck_height   )  ← axis split
#   z_step  flange outer step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,  # outer radius of the wide flange main body [m]
#     top_flange_lip_r:      float,  # outer radius of the narrower top lip [m]
#     top_flange_lip_height: float,  # axial height of the top lip [m]
#     top_flange_height:     float,  # total axial height of the flange [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids minus one bore:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Bottom face is closed (inner profile sweeps to r=0 at z=0).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin).
#          The neck inner bore (neck_inner_r at cone-axis offset) is cut through.

#       3. Upper shell [z4 → z5]
#          revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     # Inner profile descends to (0, 0) to close the bottom face.
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # inner (downward) — bore extends all the way to z=0
#         # Profile closes (neck_inner_r, 0) → (cone_bottom_outer_r, 0) = annular bottom floor
#         (neck_inner_r,        z4),
#         (neck_inner_r,        0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin; neck inner bore cut through at cone-axis offset.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         # outer (upward)
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         # closed top + solid interior: descend along axis, no diagonal
#         (0,                   z5),    # solid top face
#         (0,                   z4),    # down the axis (no surface when revolved)
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(acs)





























































# """
# Parametric SFR above-core structure (ACS) — ZLP V3.

# Three bodies unioned, one bore cut:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              bottom face CLOSED (profile goes to r=0 at z=0)

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_flange_outer_r, origin)
#                              with neck inner bore cut through at cone-axis offset

#   Upper shell   [z4 → z5]   flange (stepped outer profile, inner bore open)
#                              revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face (closed)
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = flange bottom ( = z3 + neck_height   )  ← axis split
#   z_step  flange outer step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,  # outer radius of the wide flange main body [m]
#     top_flange_lip_r:      float,  # outer radius of the narrower top lip [m]
#     top_flange_lip_height: float,  # axial height of the top lip [m]
#     top_flange_height:     float,  # total axial height of the flange [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids minus one bore:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Bottom face is closed (inner profile sweeps to r=0 at z=0).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin).
#          The neck inner bore (neck_inner_r at cone-axis offset) is cut through.

#       3. Upper shell [z4 → z5]
#          revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     # Inner profile descends to (0, 0) to close the bottom face.
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # inner (downward)
#         (neck_inner_r,        z4),
#         (neck_inner_r,        z2),
#         (cone_bottom_inner_r, z1),
#         (cone_bottom_inner_r, 0),
#         (0,                   0),    # ← closes bottom face
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin; neck inner bore cut through at cone-axis offset.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         # outer (upward)
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         # closed top + solid interior: descend along axis, no diagonal
#         (0,                   z5),    # solid top face
#         (0,                   z4),    # down the axis (no surface when revolved)
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(acs)












































# """
# Parametric SFR above-core structure (ACS) — ZLP V3.

# Three bodies unioned, one bore cut:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              bottom face CLOSED (profile goes to r=0 at z=0)

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_flange_outer_r, origin)
#                              with neck inner bore cut through at cone-axis offset

#   Upper shell   [z4 → z5]   flange (stepped outer profile, inner bore open)
#                              revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face (closed)
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = flange bottom ( = z3 + neck_height   )  ← axis split
#   z_step  flange outer step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,  # outer radius of the wide flange main body [m]
#     top_flange_lip_r:      float,  # outer radius of the narrower top lip [m]
#     top_flange_lip_height: float,  # axial height of the top lip [m]
#     top_flange_inner_r:    float,  # inner bore radius of the flange [m]
#     top_flange_height:     float,  # total axial height of the flange [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids minus one bore:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Bottom face is closed (inner profile sweeps to r=0 at z=0).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin).
#          The neck inner bore (neck_inner_r at cone-axis offset) is cut through.

#       3. Upper shell [z4 → z5]
#          Flange (stepped outer profile, inner bore = top_flange_inner_r),
#          revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_r   > top_flange_inner_r,  "top_flange_lip_r > top_flange_inner_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (top_flange_inner_r > neck_outer_r,         "top_flange_inner_r > neck_outer_r"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     # Inner profile descends to (0, 0) to close the bottom face.
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # inner (downward)
#         (neck_inner_r,        z4),
#         (neck_inner_r,        z2),
#         (cone_bottom_inner_r, z1),
#         (cone_bottom_inner_r, 0),
#         (0,                   0),    # ← closes bottom face
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin; neck inner bore cut through at cone-axis offset.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         # outer (upward)
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         # top closed: goes to axis, then conical inner surface back down
#         (0,                   z5),              # solid top face
#         (top_flange_inner_r,  z4),              # inner bore open only at z4
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_inner_r    = 1.6445,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(acs)






































# """
# Parametric SFR above-core structure (ACS) — ZLP V3.

# Three bodies unioned, one bore cut:

#   Lower shell   [0  → z4]   bottom ring + cone + collar + neck
#                              revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                              bottom face CLOSED (profile goes to r=0 at z=0)

#   Closing plate [z4-cp_h → z4]
#                              solid disk (outer r = top_flange_outer_r, origin)
#                              with neck inner bore cut through at cone-axis offset

#   Upper shell   [z4 → z5]   flange (stepped outer profile, inner bore open)
#                              revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face (closed)
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top = flange bottom ( = z3 + neck_height   )  ← axis split
#   z_step  flange outer step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_above_core_structure()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_above_core_structure(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,  # outer radius of the wide flange main body [m]
#     top_flange_lip_r:      float,  # outer radius of the narrower top lip [m]
#     top_flange_lip_height: float,  # axial height of the top lip [m]
#     top_flange_inner_r:    float,  # inner bore radius of the flange [m]
#     top_flange_height:     float,  # total axial height of the flange [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     closing_plate_height: float,   # thickness of the plate at z4 [m]

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the above-core structure as the union of three solids minus one bore:

#       1. Lower shell [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.
#          Bottom face is closed (inner profile sweeps to r=0 at z=0).

#       2. Closing plate [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin).
#          The neck inner bore (neck_inner_r at cone-axis offset) is cut through.

#       3. Upper shell [z4 → z5]
#          Flange (stepped outer profile, inner bore = top_flange_inner_r),
#          revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_r   > top_flange_inner_r,  "top_flange_lip_r > top_flange_inner_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (top_flange_inner_r > neck_outer_r,         "top_flange_inner_r > neck_outer_r"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     # Inner profile descends to (0, 0) to close the bottom face.
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         # inner (downward)
#         (neck_inner_r,        z4),
#         (neck_inner_r,        z2),
#         (cone_bottom_inner_r, z1),
#         (cone_bottom_inner_r, 0),
#         (0,                   0),    # ← closes bottom face
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk at origin; neck inner bore cut through at cone-axis offset.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         # outer (upward)
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         # inner (downward)
#         (top_flange_inner_r,  z5),
#         (top_flange_inner_r,  z4),
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     acs = create_above_core_structure(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_inner_r    = 1.6445,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(acs)


































































# """
# Parametric SFR core support barrel — ZLP V3.

# Three bodies unioned together:

#   Lower shell  [0  → z4]  bottom ring + cone + collar + neck
#                            revolved around (cone_axis_offset_x, cone_axis_offset_y)

#   Closing plate[z4-closing_plate_height → z4]
#                            solid disk (outer r = top_flange_outer_r, centered at origin)
#                            with neck bore cut through at cone-axis offset
#                            — closes the gap between the offset neck and the flange

#   Upper shell  [z4 → z5]  flange only (stepped outer profile)
#                            revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top  = flange bottom( = z3 + neck_height   )  ← axis split
#   z_step  flange step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_core_support_barrel()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_core_support_barrel(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,
#     top_flange_lip_r:      float,
#     top_flange_lip_height: float,
#     top_flange_height:     float,

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ─────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     # Solid disk centered on the upper axis that closes the gap between the
#     # offset neck top and the flange bottom at z4.
#     # Only the neck inner bore is left open (fluid passage).
#     closing_plate_height: float,

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the core support barrel as the union of three solids:

#       1. Lower shell  [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.

#       2. Closing plate  [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin) with
#          the neck inner bore cut through at the cone-axis offset position.
#          This closes the lateral gap between the offset neck and the flange.

#       3. Upper shell  [z4 → z5]
#          Flange (stepped outer profile), revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     lower_pts = [
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         (neck_inner_r,        z4),
#         (neck_inner_r,        z2),
#         (cone_bottom_inner_r, z1),
#         (cone_bottom_inner_r, 0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Solid disk centered on the upper axis (origin), no holes.
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         (0,                   z5),    # top face — solid, no bore
#         (0,                   z4),    # bottom centre — closes back to start
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes ─────────────────────────────────────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     barrel = create_core_support_barrel(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,   # adjust to match drawing
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(barrel)


































































































# """
# Parametric SFR core support barrel — ZLP V3.

# Two bodies unioned together:

#   Lower shell  [0  → z4]  bottom ring + cone + collar + neck
#                            revolved around (cone_axis_offset_x, cone_axis_offset_y)
#                            BOTTOM FACE CLOSED (profile goes to r=0 at z=0)

#   Upper shell  [z4 → z5]  flange (stepped outer profile, inner bore open)
#                            revolved around (0, 0)

# The lateral gap between the offset lower section and the flange is
# intentional and left open.

# Stacking order (bottom → top):
#   z=0     bottom face  ← closed solid floor
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top  = flange bottom( = z3 + neck_height   )  ← axis split
#   z_step  flange step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_core_support_barrel()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_core_support_barrel(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,  # outer radius of the wide flange main body [m]
#     top_flange_lip_r:      float,  # outer radius of the narrower top lip [m]
#     top_flange_lip_height: float,  # axial height of the top lip [m]
#     top_flange_inner_r:    float,  # inner bore radius of the flange [m]
#     top_flange_height:     float,  # total axial height of the flange [m]

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ─────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the core support barrel as the union of two solids:

#       Lower shell [0 → z4]
#         bottom ring + cone + collar + neck, revolved around the cone axis.
#         The bottom face is CLOSED: the inner profile runs to r=0 at z=0,
#         creating a solid floor.

#       Upper shell [z4 → z5]
#         Flange (stepped outer profile, inner bore = top_flange_inner_r),
#         revolved around the upper axis (0, 0).
#         The lateral gap between the offset lower section and the flange
#         is intentionally left open.
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_r   > top_flange_inner_r,  "top_flange_lip_r > top_flange_inner_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (top_flange_inner_r > neck_outer_r,         "top_flange_inner_r > neck_outer_r"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell — bottom ring + cone + collar + neck ───────────────
#     # The inner profile descends from neck_inner_r at z4 down to
#     # cone_bottom_inner_r at z1, then to (0, 0) — closing the bottom face.
#     lower_pts = [
#         # outer surface (upward)
#         (cone_bottom_outer_r, 0),    # A  bottom outer edge
#         (cone_bottom_outer_r, z1),   # B  bottom ring outer top
#         (neck_outer_r,        z2),   # C  cone outer top = collar outer bottom
#         (collar_outer_r,      z2),   # D  collar outer bottom
#         (collar_outer_r,      z3),   # E  collar outer top = neck outer bottom
#         (neck_outer_r,        z3),   # F  neck outer bottom
#         (neck_outer_r,        z4),   # G  neck outer top
#         # inner surface (downward)
#         (neck_inner_r,        z4),   # H  neck inner top
#         (neck_inner_r,        z2),   # I  neck + collar inner bore
#         (cone_bottom_inner_r, z1),   # J  cone inner bottom
#         (cone_bottom_inner_r, 0),    # K  inner bottom edge
#         (0,                   0),    # L  axis at z=0 → closes bottom face
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Upper shell — flange ───────────────────────────────────────────
#     upper_pts = [
#         # outer surface (upward)
#         (top_flange_outer_r,  z4),      # a  flange outer bottom
#         (top_flange_outer_r,  z_step),  # b  flange main body top
#         (top_flange_lip_r,    z_step),  # c  step inward to lip
#         (top_flange_lip_r,    z5),      # d  flange top outer
#         # inner surface (downward)
#         (top_flange_inner_r,  z5),      # e  flange top inner
#         (top_flange_inner_r,  z4),      # f  flange bottom inner → closes to a
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 3. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(lower_solid)

#     # ── 4. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 5. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     barrel = create_core_support_barrel(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_inner_r    = 1.6445,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(barrel)


































































# """
# Parametric SFR core support barrel — ZLP V3.

# Three bodies unioned together:

#   Lower shell  [0  → z4]  bottom ring + cone + collar + neck
#                            revolved around (cone_axis_offset_x, cone_axis_offset_y)

#   Closing plate[z4-closing_plate_height → z4]
#                            solid disk (outer r = top_flange_outer_r, centered at origin)
#                            with neck bore cut through at cone-axis offset
#                            — closes the gap between the offset neck and the flange

#   Upper shell  [z4 → z5]  flange only (stepped outer profile)
#                            revolved around (0, 0)

# Stacking order (bottom → top):
#   z=0     bottom face
#   z1      bottom ring top          ( = bottom_ring_height )
#   z2      cone top = collar bottom ( = z1 + cone_height   )
#   z3      collar top = neck bottom ( = z2 + collar_height )
#   z4      neck top  = flange bottom( = z3 + neck_height   )  ← axis split
#   z_step  flange step
#   z5      flange top

# All geometry dimensions are required — no defaults.
# Single public function:  create_core_support_barrel()
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_core_support_barrel(
#     # ── Flange — upper axis (0, 0) ────────────────────────────────────────
#     top_flange_outer_r:    float,
#     top_flange_lip_r:      float,
#     top_flange_lip_height: float,
#     top_flange_inner_r:    float,
#     top_flange_height:     float,

#     # ── Neck — lower / cone axis ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar — lower / cone axis ─────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness (uniform) ──────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring — lower / cone axis ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Closing plate ─────────────────────────────────────────────────────
#     # Solid disk centered on the upper axis that closes the gap between the
#     # offset neck top and the flange bottom at z4.
#     # Only the neck inner bore is left open (fluid passage).
#     closing_plate_height: float,

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the core support barrel as the union of three solids:

#       1. Lower shell  [0 → z4]
#          bottom ring + cone + collar + neck, revolved around the cone axis.

#       2. Closing plate  [z4 - closing_plate_height → z4]
#          Solid disk (outer r = top_flange_outer_r, centered at origin) with
#          the neck inner bore cut through at the cone-axis offset position.
#          This closes the lateral gap between the offset neck and the flange.

#       3. Upper shell  [z4 → z5]
#          Flange (stepped outer profile), revolved around the upper axis (0, 0).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,    "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_r   > top_flange_inner_r,  "top_flange_lip_r > top_flange_inner_r"),
#         (top_flange_lip_height < top_flange_height, "top_flange_lip_height < top_flange_height"),
#         (top_flange_inner_r > neck_outer_r,         "top_flange_inner_r > neck_outer_r"),
#         (collar_outer_r >= neck_outer_r,            "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                     "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,        "cone_bottom_outer_r > neck_outer_r"),
#         (closing_plate_height > 0,                  "closing_plate_height > 0"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height
#     z3     = z2 + collar_height
#     z4     = z3 + neck_height                                   # ← axis split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell ────────────────────────────────────────────────────
#     lower_pts = [
#         (cone_bottom_outer_r, 0),
#         (cone_bottom_outer_r, z1),
#         (neck_outer_r,        z2),
#         (collar_outer_r,      z2),
#         (collar_outer_r,      z3),
#         (neck_outer_r,        z3),
#         (neck_outer_r,        z4),
#         (neck_inner_r,        z4),
#         (neck_inner_r,        z2),
#         (cone_bottom_inner_r, z1),
#         (cone_bottom_inner_r, 0),
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Closing plate ──────────────────────────────────────────────────
#     # Full disk centered on the upper axis (origin), sitting from
#     # z4 - closing_plate_height to z4.  The neck inner bore is the only hole.
#     EPS = 1e-4
#     closing_plate = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height)
#         .circle(top_flange_outer_r)
#         .extrude(closing_plate_height)
#     )
#     neck_bore_cutter = (
#         cq.Workplane("XY")
#         .workplane(offset=z4 - closing_plate_height - EPS)
#         .circle(neck_inner_r)
#         .extrude(closing_plate_height + 2 * EPS)
#         .translate((cone_axis_offset_x, cone_axis_offset_y, 0))
#     )
#     closing_plate = closing_plate.cut(neck_bore_cutter)

#     # ── 3. Upper shell: flange ────────────────────────────────────────────
#     upper_pts = [
#         (top_flange_outer_r,  z4),
#         (top_flange_outer_r,  z_step),
#         (top_flange_lip_r,    z_step),
#         (top_flange_lip_r,    z5),
#         (top_flange_inner_r,  z5),
#         (top_flange_inner_r,  z4),
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 4. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(closing_plate).union(lower_solid)

#     # ── 5. Flow holes ─────────────────────────────────────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 6. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     barrel = create_core_support_barrel(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_inner_r    = 1.6445,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         closing_plate_height  = 0.050,   # adjust to match drawing
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(barrel)

















    











































# """
# Parametric SFR core support barrel — ZLP V3.

# Upper shell [z4 → z5] — FLANGE ONLY, revolved around (0, 0)
# Lower shell [0  → z4] — bottom ring + cone + collar + neck,
#                          revolved around (cone_axis_offset_x, cone_axis_offset_y)

# Stacking order (bottom → top):
#   z=0    bottom face
#   z1     bottom ring top          ( = bottom_ring_height )
#   z2     cone top = collar bottom ( = z1 + cone_height   )
#   z3     collar top = neck bottom ( = z2 + collar_height )
#   z4     neck top  = flange bottom( = z3 + neck_height   )  ← axis split
#   z_step flange step              ( = z4 + flange_h - lip_h )
#   z5     flange top               ( = z4 + flange_height )

# All geometry dimensions are required — no defaults.
# Single public function:  create_core_support_barrel()

# Key ESFR-SMART nominal dimensions (metres):
#     top_flange_outer_r    = 1.843    (Ø3686 / 2)
#     top_flange_lip_r      = 1.7055   (Ø3411 / 2)
#     top_flange_lip_height = 0.250
#     top_flange_inner_r    = 1.6445   (Ø3289 / 2)
#     top_flange_height     = 1.008
#     neck_outer_r          = 1.1085   (Ø2217 / 2)
#     neck_height           = 0.569
#     collar_outer_r        = 1.1085
#     collar_height         = 0.092
#     wall_t                = 0.025
#     cone_height           = 2.429
#     cone_bottom_outer_r   = 1.403    (Ø2806 / 2)
#     bottom_ring_height    = 0.498
#     cone_axis_offset_x    = 0.6056   (605.6 mm)
#     cone_axis_offset_y    = 0.0
# """

# from __future__ import annotations
# import math
# import cadquery as cq

# from profile_from_straight_connections import create_profile_from_straight_connections
# from utils import revolve_profile


# def _dedup(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
#     EPS = 1e-9
#     out = [pts[0]]
#     for p in pts[1:]:
#         if abs(p[0] - out[-1][0]) > EPS or abs(p[1] - out[-1][1]) > EPS:
#             out.append(p)
#     return out


# def _revolve_closed(pts: list[tuple[float, float]]) -> cq.Workplane:
#     profile = create_profile_from_straight_connections(_dedup(pts), plane="XZ", closed=True)
#     return revolve_profile(profile, angle=360, axis="Z")


# def create_core_support_barrel(
#     # ── Flange (upper axis) ───────────────────────────────────────────────
#     top_flange_outer_r:    float,
#     top_flange_lip_r:      float,
#     top_flange_lip_height: float,
#     top_flange_inner_r:    float,
#     top_flange_height:     float,

#     # ── Neck (lower / cone axis) ──────────────────────────────────────────
#     neck_outer_r: float,
#     neck_height:  float,

#     # ── Collar (lower / cone axis) ────────────────────────────────────────
#     collar_outer_r: float,
#     collar_height:  float,

#     # ── Wall thickness ────────────────────────────────────────────────────
#     wall_t: float,

#     # ── Cone + bottom ring (lower / cone axis) ────────────────────────────
#     cone_bottom_outer_r: float,
#     cone_height:         float,
#     bottom_ring_height:  float,

#     # ── Cone axis offset ──────────────────────────────────────────────────
#     cone_axis_offset_x: float,
#     cone_axis_offset_y: float,

#     # ── Optional flow holes on the cone ───────────────────────────────────
#     flow_hole_groups: list[dict] | None = None,

#     # ── Global position ───────────────────────────────────────────────────
#     z_bottom: float = 0.0,

# ) -> cq.Workplane:
#     """
#     Build the core support barrel as the union of:

#       Lower shell [0 → z4]  bottom ring + cone + collar + neck
#                              revolved around Z at (cone_axis_offset_x, cone_axis_offset_y)

#       Upper shell [z4 → z5] flange only (stepped outer profile)
#                              revolved around Z at (0, 0)

#     The two shells share the z4 plane (neck top / flange bottom).
#     """
#     # ── Validate ─────────────────────────────────────────────────────────
#     checks = [
#         (top_flange_outer_r > top_flange_lip_r,   "top_flange_outer_r > top_flange_lip_r"),
#         (top_flange_lip_r   > top_flange_inner_r, "top_flange_lip_r > top_flange_inner_r"),
#         (top_flange_lip_height < top_flange_height,"top_flange_lip_height < top_flange_height"),
#         (top_flange_inner_r > neck_outer_r,        "top_flange_inner_r > neck_outer_r"),
#         (collar_outer_r >= neck_outer_r,           "collar_outer_r >= neck_outer_r"),
#         (neck_outer_r > wall_t,                    "neck_outer_r > wall_t"),
#         (cone_bottom_outer_r > neck_outer_r,       "cone_bottom_outer_r > neck_outer_r"),
#     ]
#     for ok, msg in checks:
#         if not ok:
#             raise ValueError(f"Validation failed: {msg}")
#     for name, val in [
#         ("top_flange_height", top_flange_height), ("top_flange_lip_height", top_flange_lip_height),
#         ("neck_height", neck_height), ("collar_height", collar_height),
#         ("cone_height", cone_height), ("bottom_ring_height", bottom_ring_height),
#         ("wall_t", wall_t),
#     ]:
#         if val <= 0:
#             raise ValueError(f"{name} must be > 0")

#     # ── Z levels ──────────────────────────────────────────────────────────
#     z1     = bottom_ring_height
#     z2     = z1 + cone_height                                   # cone top = collar bottom
#     z3     = z2 + collar_height                                 # collar top = neck bottom
#     z4     = z3 + neck_height                                   # neck top = flange bottom  ← split
#     z_step = z4 + (top_flange_height - top_flange_lip_height)
#     z5     = z4 + top_flange_height

#     neck_inner_r        = neck_outer_r - wall_t
#     cone_bottom_inner_r = cone_bottom_outer_r - wall_t

#     # ── 1. Lower shell: bottom ring + cone + collar + neck ────────────────
#     lower_pts = [
#         # outer (upward)
#         (cone_bottom_outer_r, 0),    # A  bottom outer
#         (cone_bottom_outer_r, z1),   # B  bottom ring outer top
#         (neck_outer_r,        z2),   # C  cone outer top = collar outer bottom
#         (collar_outer_r,      z2),   # D  collar outer bottom (no-op if == C)
#         (collar_outer_r,      z3),   # E  collar outer top = neck outer bottom
#         (neck_outer_r,        z3),   # F  neck outer bottom (no-op if == E)
#         (neck_outer_r,        z4),   # G  neck outer top
#         # inner (downward)
#         (neck_inner_r,        z4),   # H  neck inner top
#         (neck_inner_r,        z2),   # I  neck + collar inner (straight bore)
#         (cone_bottom_inner_r, z1),   # J  cone inner bottom
#         (cone_bottom_inner_r, 0),    # K  bottom inner  → closes to A
#     ]
#     lower_solid = _revolve_closed(lower_pts)
#     if cone_axis_offset_x != 0.0 or cone_axis_offset_y != 0.0:
#         lower_solid = lower_solid.translate((cone_axis_offset_x, cone_axis_offset_y, 0))

#     # ── 2. Upper shell: flange only ───────────────────────────────────────
#     upper_pts = [
#         # outer (upward)
#         (top_flange_outer_r,  z4),      # a  flange outer bottom
#         (top_flange_outer_r,  z_step),  # b  flange main body top
#         (top_flange_lip_r,    z_step),  # c  step inward to lip
#         (top_flange_lip_r,    z5),      # d  flange top outer
#         # inner (downward)
#         (top_flange_inner_r,  z5),      # e  flange top inner
#         (top_flange_inner_r,  z4),      # f  flange bottom inner  → closes to a
#     ]
#     upper_solid = _revolve_closed(upper_pts)

#     # ── 3. Union ──────────────────────────────────────────────────────────
#     solid = upper_solid.union(lower_solid)

#     # ── 4. Flow holes (positions in cone-axis frame) ──────────────────────
#     cx, cy = cone_axis_offset_x, cone_axis_offset_y
#     if flow_hole_groups:
#         for group in flow_hole_groups:
#             hole_r  = float(group["hole_r"])
#             z_c     = float(group["z_center"])
#             n       = int(group["n_holes"])
#             start_a = float(group.get("start_angle_deg", 0.0))
#             if z_c < 0 or z_c > z5:
#                 raise ValueError(f"z_center={z_c} outside [0, {z5:.4f}]")
#             if z_c <= z1:
#                 r_mid = cone_bottom_outer_r - wall_t / 2.0
#             elif z_c <= z2:
#                 frac  = (z_c - z1) / (z2 - z1)
#                 r_out = cone_bottom_outer_r + frac * (neck_outer_r - cone_bottom_outer_r)
#                 r_mid = r_out - wall_t / 2.0
#             else:
#                 r_mid = neck_outer_r - wall_t / 2.0
#             for i in range(n):
#                 a  = math.radians(start_a + 360.0 * i / n)
#                 hx = cx + r_mid * math.cos(a)
#                 hy = cy + r_mid * math.sin(a)
#                 cutter = (
#                     cq.Workplane("XY").workplane(offset=z_c)
#                     .circle(hole_r).extrude(wall_t * 2)
#                     .translate((hx, hy, -wall_t))
#                 )
#                 solid = solid.cut(cutter)

#     # ── 5. Final z translation ────────────────────────────────────────────
#     if z_bottom != 0.0:
#         solid = solid.translate((0, 0, z_bottom))

#     return solid


# if __name__ == "__main__":
#     from ocp_vscode import show

#     barrel = create_core_support_barrel(
#         top_flange_outer_r    = 1.843,
#         top_flange_lip_r      = 1.7055,
#         top_flange_lip_height = 0.250,
#         top_flange_inner_r    = 1.6445,
#         top_flange_height     = 1.008,
#         neck_outer_r          = 1.1085,
#         neck_height           = 0.569,
#         collar_outer_r        = 1.1085,
#         collar_height         = 0.092,
#         wall_t                = 0.025,
#         cone_height           = 2.429,
#         cone_bottom_outer_r   = 1.403,
#         bottom_ring_height    = 0.498,
#         cone_axis_offset_x    = 0.6056,
#         cone_axis_offset_y    = 0.0,
#     )
#     show(barrel)



