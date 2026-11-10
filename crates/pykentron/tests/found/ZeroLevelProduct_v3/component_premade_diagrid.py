"""
Parametric diagrid — final version.

Hollow short cylinder with thin walls (default 30 mm on all sides) closed
top and bottom, with lateral bosses for pump-elbow connection. The bores
pierce the lateral wall and open into the interior cavity.

Key changes vs early versions:
  - `wall_t` split into `wall_t_side`, `wall_t_top`, `wall_t_bottom` so the
    plates can be thinner than the lateral wall (or vice versa) without
    blocking the bore.
  - Bore length is computed analytically so the LATERAL edges of the
    bore — not just the centerline — pierce the inner wall, avoiding the
    "vertical strip" rendering artifact.
  - Bosses are inset into the outer wall before being extruded radially
    outward, so the union with the curved wall is gap-free.
  - Validation: raises a clear error if the top/bottom plate would clip
    the bore, or if the boss vertically extends past the disc faces.
"""

from __future__ import annotations
import math
import cadquery as cq


def create_diagrid(
    diameter:               float,
    thickness:              float,
    z_bottom:               float = 0.0,
    wall_t:                 float | None = None,   # legacy alias for all three
    wall_t_side:            float = 0.030,
    wall_t_top:             float = 0.030,
    wall_t_bottom:          float = 0.030,
    open_top:               bool  = False,
    open_bottom:            bool  = False,
    nozzle_boss_angles_deg: list[float] | None = None,
    nozzle_z_abs:           float | None = None,
    nozzle_r_bore:          float = 0.230,
    nozzle_r_boss:          float = 0.301,
    nozzle_boss_height:     float = 0.0775,
) -> cq.Workplane:

    # Legacy support
    if wall_t is not None:
        wall_t_side = wall_t_top = wall_t_bottom = wall_t

    if diameter <= 0 or thickness <= 0:
        raise ValueError("diameter and thickness must be > 0")
    if wall_t_side <= 0:
        raise ValueError("wall_t_side must be > 0")
    if wall_t_top < 0 or wall_t_bottom < 0:
        raise ValueError("wall_t_top and wall_t_bottom must be ≥ 0")

    radius_outer = diameter / 2.0
    radius_inner = radius_outer - wall_t_side
    if radius_inner <= 0:
        raise ValueError(
            f"wall_t_side ({wall_t_side}) is too large for diameter ({diameter})."
        )

    z_top = z_bottom + thickness

    cavity_z_bottom = z_bottom if open_bottom else z_bottom + wall_t_bottom
    cavity_z_top    = z_top    if open_top    else z_top    - wall_t_top
    cavity_height   = cavity_z_top - cavity_z_bottom
    if cavity_height <= 0:
        raise ValueError(
            f"wall_t_top + wall_t_bottom ({wall_t_top + wall_t_bottom}) "
            f"is too large for thickness ({thickness})."
        )

    # Bore obstruction check
    if nozzle_z_abs is not None and nozzle_boss_angles_deg:
        bore_z_low  = nozzle_z_abs - nozzle_r_bore
        bore_z_high = nozzle_z_abs + nozzle_r_bore
        if (not open_top) and bore_z_high > cavity_z_top:
            raise ValueError(
                f"Top plate (Z = [{cavity_z_top:.4f}, {z_top:.4f}]) "
                f"overlaps the bore (Z = [{bore_z_low:.4f}, "
                f"{bore_z_high:.4f}]). Reduce wall_t_top to ≤ "
                f"{(z_top - bore_z_high):.4f} m or move the bore down."
            )
        if (not open_bottom) and bore_z_low < cavity_z_bottom:
            raise ValueError(
                f"Bottom plate (Z = [{z_bottom:.4f}, {cavity_z_bottom:.4f}]) "
                f"overlaps the bore (Z = [{bore_z_low:.4f}, "
                f"{bore_z_high:.4f}]). Reduce wall_t_bottom to ≤ "
                f"{(bore_z_low - z_bottom):.4f} m or move the bore up."
            )

    # Build hollow body
    disc = (cq.Workplane("XY")
            .workplane(offset=z_bottom)
            .circle(radius_outer)
            .extrude(thickness))

    cavity = (cq.Workplane("XY")
              .workplane(offset=cavity_z_bottom)
              .circle(radius_inner)
              .extrude(cavity_height))

    disc = disc.cut(cavity)

    # Bosses + bores
    if nozzle_boss_angles_deg and nozzle_z_abs is not None:

        if nozzle_r_boss >= radius_outer:
            raise ValueError(
                f"nozzle_r_boss ({nozzle_r_boss}) must be smaller than "
                f"diagrid outer radius ({radius_outer})."
            )
        if nozzle_r_bore >= radius_inner:
            raise ValueError(
                f"nozzle_r_bore ({nozzle_r_bore}) must be smaller than "
                f"diagrid inner radius ({radius_inner})."
            )

        inset_min  = radius_outer - math.sqrt(radius_outer**2 - nozzle_r_boss**2)
        inset_safe = inset_min + 0.005
        boss_total_height = inset_safe + nozzle_boss_height

        margin = nozzle_r_boss
        if not (z_bottom + margin < nozzle_z_abs < z_top - margin):
            raise ValueError(
                f"nozzle_z_abs={nozzle_z_abs:.4f} m clips the diagrid face. "
                f"Must be in ({z_bottom + margin:.4f}, {z_top - margin:.4f})."
            )

        # Bore length so the lateral edges of the bore also pierce the inner wall.
        L_min = (
            (radius_outer + nozzle_boss_height)
            - math.sqrt(radius_inner**2 - nozzle_r_bore**2)
        )
        bore_length = L_min + 0.010

        for theta in nozzle_boss_angles_deg:
            rad = math.radians(theta)

            outward = cq.Vector( math.cos(rad),  math.sin(rad), 0.0)
            inward  = cq.Vector(-math.cos(rad), -math.sin(rad), 0.0)
            tangent = cq.Vector(-math.sin(rad),  math.cos(rad), 0.0)

            base_x = (radius_outer - inset_safe) * math.cos(rad)
            base_y = (radius_outer - inset_safe) * math.sin(rad)
            base_origin = cq.Vector(base_x, base_y, nozzle_z_abs)

            boss_plane = cq.Plane(
                origin=base_origin, xDir=tangent, normal=outward
            )
            boss = (cq.Workplane(boss_plane)
                    .circle(nozzle_r_boss)
                    .extrude(boss_total_height))
            disc = disc.union(boss)

            outer_face_origin = cq.Vector(
                (radius_outer + nozzle_boss_height) * math.cos(rad),
                (radius_outer + nozzle_boss_height) * math.sin(rad),
                nozzle_z_abs,
            )
            bore_plane = cq.Plane(
                origin=outer_face_origin, xDir=tangent, normal=inward
            )
            bore = (cq.Workplane(bore_plane)
                    .circle(nozzle_r_bore)
                    .extrude(bore_length))
            disc = disc.cut(bore)

    return disc.clean()


if __name__ == "__main__":
    from ocp_vscode import show

    diagrid = create_diagrid(
        diameter   = 4.660,
        thickness  = 1.050,
        z_bottom   = -0.460,
    )
    show(diagrid)







































# FUNCIONA BIEN

# """
# Parametric diagrid — ZLP V11.

# Fix vs V10:
# ─────────────
# 7) Top/bottom plates no longer overlap the boss bores.

#    V10 used a single `wall_t` for the lateral wall AND for the top and
#    bottom plates. With wall_t = 0.300 m and the boss centered at
#    Z ≈ 0.247 m (radius 0.230 m), the top plate occupied Z = [0.290,
#    0.590] which OVERLAPPED the upper part of the bore (which reached
#    up to Z = 0.477). Result: ≈ 187 mm of the bore was visually open
#    on the outer face but BLOCKED by the top plate behind it — no
#    fluid passage.

#    Fix: split into `wall_t_side` (lateral) and `wall_t_top` /
#    `wall_t_bottom` (plates). The lateral wall remains thick (300 mm)
#    for structural strength; the top and bottom plates are thinner
#    (default 50 mm), and we additionally REQUIRE that the plates do
#    NOT intersect the boss bores when the boss positions are known.
#    We raise an error if they would.

#    For backward compatibility we keep `wall_t` as an alias: if given,
#    it sets all three (side, top, bottom) to the same value.
# """

# from __future__ import annotations
# import math
# import cadquery as cq


# def create_diagrid(
#     diameter:               float,
#     thickness:              float,
#     z_bottom:               float = 0.0,
#     wall_t:                 float | None = None,   # legacy alias
#     wall_t_side:            float = 0.030,
#     wall_t_top:             float = 0.030,
#     wall_t_bottom:          float = 0.030,
#     open_top:               bool  = False,
#     open_bottom:            bool  = False,
#     nozzle_boss_angles_deg: list[float] | None = None,
#     nozzle_z_abs:           float | None = None,
#     nozzle_r_bore:          float = 0.230,
#     nozzle_r_boss:          float = 0.301,
#     nozzle_boss_height:     float = 0.0775,
# ) -> cq.Workplane:

#     # Legacy support: a single wall_t overrides all three.
#     if wall_t is not None:
#         wall_t_side = wall_t_top = wall_t_bottom = wall_t

#     if diameter <= 0 or thickness <= 0:
#         raise ValueError("diameter and thickness must be > 0")
#     if wall_t_side <= 0:
#         raise ValueError("wall_t_side must be > 0")
#     if wall_t_top < 0 or wall_t_bottom < 0:
#         raise ValueError("wall_t_top and wall_t_bottom must be ≥ 0")

#     radius_outer = diameter / 2.0
#     radius_inner = radius_outer - wall_t_side
#     if radius_inner <= 0:
#         raise ValueError(
#             f"wall_t_side ({wall_t_side}) is too large: inner cavity "
#             f"radius would be {radius_inner}."
#         )

#     z_top = z_bottom + thickness

#     # ── Cavity vertical extent ────────────────────────────────────────
#     cavity_z_bottom = z_bottom if open_bottom else z_bottom + wall_t_bottom
#     cavity_z_top    = z_top    if open_top    else z_top    - wall_t_top
#     cavity_height   = cavity_z_top - cavity_z_bottom
#     if cavity_height <= 0:
#         raise ValueError(
#             f"wall_t_top + wall_t_bottom ({wall_t_top + wall_t_bottom}) "
#             f"is too large for thickness ({thickness}): top and bottom "
#             "plates would overlap."
#         )

#     # ── Bore obstruction check ────────────────────────────────────────
#     # If we know where the bores will be (nozzle_z_abs known), make sure
#     # the top and bottom plates are NOT in the bore's Z range. Bore
#     # spans Z = [nozzle_z_abs − r_bore, nozzle_z_abs + r_bore]. The
#     # plates span Z = [z_bottom, cavity_z_bottom] and [cavity_z_top, z_top].
#     if nozzle_z_abs is not None and nozzle_boss_angles_deg:
#         bore_z_low  = nozzle_z_abs - nozzle_r_bore
#         bore_z_high = nozzle_z_abs + nozzle_r_bore
#         if (not open_top) and bore_z_high > cavity_z_top:
#             raise ValueError(
#                 f"Top plate (Z = [{cavity_z_top:.4f}, {z_top:.4f}]) "
#                 f"overlaps the bore (Z = [{bore_z_low:.4f}, "
#                 f"{bore_z_high:.4f}]). Reduce wall_t_top to ≤ "
#                 f"{(z_top - bore_z_high):.4f} m or move the bore down."
#             )
#         if (not open_bottom) and bore_z_low < cavity_z_bottom:
#             raise ValueError(
#                 f"Bottom plate (Z = [{z_bottom:.4f}, {cavity_z_bottom:.4f}]) "
#                 f"overlaps the bore (Z = [{bore_z_low:.4f}, "
#                 f"{bore_z_high:.4f}]). Reduce wall_t_bottom to ≤ "
#                 f"{(bore_z_low - z_bottom):.4f} m or move the bore up."
#             )

#     # ── Build hollow diagrid ──────────────────────────────────────────
#     disc = (cq.Workplane("XY")
#             .workplane(offset=z_bottom)
#             .circle(radius_outer)
#             .extrude(thickness))

#     cavity = (cq.Workplane("XY")
#               .workplane(offset=cavity_z_bottom)
#               .circle(radius_inner)
#               .extrude(cavity_height))

#     disc = disc.cut(cavity)

#     # ── Bosses + bores ────────────────────────────────────────────────
#     if nozzle_boss_angles_deg and nozzle_z_abs is not None:

#         if nozzle_r_boss >= radius_outer:
#             raise ValueError(
#                 f"nozzle_r_boss ({nozzle_r_boss}) must be smaller than "
#                 f"diagrid outer radius ({radius_outer})."
#             )
#         if nozzle_r_bore >= radius_inner:
#             raise ValueError(
#                 f"nozzle_r_bore ({nozzle_r_bore}) must be smaller than "
#                 f"diagrid inner radius ({radius_inner})."
#             )

#         inset_min  = radius_outer - math.sqrt(radius_outer**2 - nozzle_r_boss**2)
#         inset_safe = inset_min + 0.005
#         boss_total_height = inset_safe + nozzle_boss_height

#         # Boss vertical clearance check (against full disc, not cavity).
#         margin = nozzle_r_boss
#         if not (z_bottom + margin < nozzle_z_abs < z_top - margin):
#             raise ValueError(
#                 f"nozzle_z_abs={nozzle_z_abs:.4f} m clips the diagrid face. "
#                 f"Must be in ({z_bottom + margin:.4f}, {z_top - margin:.4f})."
#             )

#         # Bore length so the LATERAL edges also pierce the inner wall.
#         L_min = (
#             (radius_outer + nozzle_boss_height)
#             - math.sqrt(radius_inner**2 - nozzle_r_bore**2)
#         )
#         bore_length = L_min + 0.010

#         for theta in nozzle_boss_angles_deg:
#             rad = math.radians(theta)

#             outward = cq.Vector( math.cos(rad),  math.sin(rad), 0.0)
#             inward  = cq.Vector(-math.cos(rad), -math.sin(rad), 0.0)
#             tangent = cq.Vector(-math.sin(rad),  math.cos(rad), 0.0)

#             base_x = (radius_outer - inset_safe) * math.cos(rad)
#             base_y = (radius_outer - inset_safe) * math.sin(rad)
#             base_origin = cq.Vector(base_x, base_y, nozzle_z_abs)

#             boss_plane = cq.Plane(
#                 origin=base_origin, xDir=tangent, normal=outward
#             )
#             boss = (cq.Workplane(boss_plane)
#                     .circle(nozzle_r_boss)
#                     .extrude(boss_total_height))
#             disc = disc.union(boss)

#             outer_face_origin = cq.Vector(
#                 (radius_outer + nozzle_boss_height) * math.cos(rad),
#                 (radius_outer + nozzle_boss_height) * math.sin(rad),
#                 nozzle_z_abs,
#             )
#             bore_plane = cq.Plane(
#                 origin=outer_face_origin, xDir=tangent, normal=inward
#             )
#             bore = (cq.Workplane(bore_plane)
#                     .circle(nozzle_r_bore)
#                     .extrude(bore_length))
#             disc = disc.cut(bore)

#     return disc.clean()


# # ──────────────────────────────────────────────────────────────────────
# # Demo
# # ──────────────────────────────────────────────────────────────────────
# if __name__ == "__main__":
#     from ocp_vscode import show

#     _RPV_INNER_D = 8.91;  _RPV_WALL_T = 0.05
#     _TORI_Rc = 5.245;     _TORI_rk = 0.379
#     _od = _RPV_INNER_D + 2 * _RPV_WALL_T
#     _r  = _od / 2
#     _xk = _r - _TORI_rk
#     _zc = math.sqrt((_TORI_Rc - _TORI_rk)**2 - _xk**2)
#     _HEAD_BOTTOM_Z    = _zc - _TORI_Rc
#     _PUMP_Z_BOTTOM    = _HEAD_BOTTOM_Z + 2.562
#     _DIAGRID_Z_BOTTOM = -1.702 + 1.242
#     _pump_r           = 3.369

#     arc_rad   = math.radians(105.0)
#     L_leg     = 0.600; L_inlet = 0.050; R_bend = 0.460
#     barrel_r  = 1.350 / 2
#     overshoot = 0.040 * 1.5
#     ex  = R_bend * (1.0 - math.cos(arc_rad)) + L_leg * math.sin(arc_rad)
#     ey  = L_inlet + R_bend * math.sin(arc_rad) + L_leg * math.cos(arc_rad)
#     lx  = ey + (barrel_r - overshoot)
#     ly  = -ex
#     phi = math.degrees(math.atan2(lx, _pump_r + ly))

#     boss_angles = []
#     for a in [60.0, 180.0, 300.0]:
#         boss_angles.append(a - phi)
#         boss_angles.append(a + phi)

#     # Bore Z range with nozzle_z_abs = 0.247, r_bore = 0.230:
#     #   bore spans Z = [0.017, 0.477]
#     #   wall_t_top max so it doesn't intersect bore: z_top − 0.477 = 0.113
#     #   wall_t_bottom max: 0.017 − z_bottom = 0.477
#     # We use 30 mm thin walls on all sides.
#     diagrid = create_diagrid(
#         diameter               = 4.660,
#         thickness              = 1.050,
#         z_bottom               = _DIAGRID_Z_BOTTOM,
#         wall_t_side            = 0.030,
#         wall_t_top             = 0.030,
#         wall_t_bottom          = 0.030,
#         open_top               = False,
#         open_bottom            = False,
#         nozzle_boss_angles_deg = boss_angles,
#         nozzle_z_abs           = _PUMP_Z_BOTTOM + 0.350,
#         nozzle_r_bore          = 0.230,
#         nozzle_r_boss          = 0.301,
#         nozzle_boss_height     = 0.0775,
#     )
#     show(diagrid)








































































# """
# Parametric diagrid — ZLP V7.

# Fixes vs V6:
# ─────────────
# 1) Bosses no longer "float" off the disc.
#    In V6, each boss was built on a plane *tangent* to the diagrid's outer
#    cylindrical wall (its base circle sat on the tangent line), so the
#    curved cylindrical wall only met the boss along a single vertical line
#    and left visible gaps at the sides of each boss (see screenshot).

#    Fix: build the boss starting from a plane that is INSET into the disc
#    by `_inset` metres, and extrude radially OUTWARD a longer distance
#    `boss_total_height = _inset + nozzle_boss_height`. The result is a
#    cylinder that physically pierces the curved wall and emerges as a
#    clean circular boss on the outside. The `_inset` is computed from
#    geometry so the base circle of the boss is fully contained inside the
#    diagrid wall:
#        inset_min = R - sqrt(R² - r_boss²)
#    plus a small safety margin.

# 2) Boss outer face is positioned to match the pump nozzle mouth radius.
#    The j-bent pump nozzle terminates at an absolute radius of ≈ 2.4075 m
#    (computed from elbow geometry). The boss outer face must land exactly
#    there so the two pieces meet flush. With diagrid radius = 2.330,
#    the boss needs to protrude (2.4075 − 2.330) = 0.0775 m.

# The `nozzle_boss_height` argument is now interpreted as the *protrusion*
# beyond the diagrid outer wall (NOT the total extrusion height). The inset
# is handled internally and not user-tunable.
# """

# from __future__ import annotations
# import math
# import cadquery as cq


# def create_diagrid(
#     diameter:               float,
#     thickness:              float,
#     z_bottom:               float = 0.0,
#     nozzle_boss_angles_deg: list[float] | None = None,
#     nozzle_z_abs:           float | None = None,
#     nozzle_r_bore:          float = 0.230,
#     nozzle_depth:           float = 0.300,
#     nozzle_r_boss:          float = 0.301,
#     nozzle_boss_height:     float = 0.0775,   # PROTRUSION beyond outer wall
# ) -> cq.Workplane:

#     if diameter <= 0 or thickness <= 0:
#         raise ValueError("diameter and thickness must be > 0")

#     radius = diameter / 2.0
#     z_top  = z_bottom + thickness

#     disc = (cq.Workplane("XY")
#             .workplane(offset=z_bottom)
#             .circle(radius)
#             .extrude(thickness))

#     if nozzle_boss_angles_deg and nozzle_z_abs is not None:

#         # ── Compute geometric inset so the boss base circle fits inside the
#         #    cylindrical wall (no gaps at the boss/wall fillet).
#         #
#         #    The disc is a cylinder of radius R. A cylinder of radius r_boss
#         #    aimed radially is fully contained in the disc wall (in the
#         #    sense that its base circle lies inside the wall) iff its base
#         #    plane is at radial distance ≤ sqrt(R² − r_boss²) from the disc
#         #    axis. So the minimum inset of the base plane from the outer
#         #    wall is:
#         #        inset_min = R − sqrt(R² − r_boss²)
#         #
#         #    We add a small safety margin so the union with the disc
#         #    produces a clean, watertight body.
#         if nozzle_r_boss >= radius:
#             raise ValueError(
#                 f"nozzle_r_boss ({nozzle_r_boss}) must be smaller than "
#                 f"diagrid radius ({radius})."
#             )
#         inset_min  = radius - math.sqrt(radius**2 - nozzle_r_boss**2)
#         inset_safe = inset_min + 0.005          # 5 mm safety margin

#         # Total radial extrusion length of the boss: from the inset plane
#         # outward, past the wall, to a protrusion of `nozzle_boss_height`
#         # beyond the outer wall.
#         boss_total_height = inset_safe + nozzle_boss_height

#         # Vertical clearance check: the boss must not clip the top/bottom
#         # faces of the diagrid.
#         margin = nozzle_r_boss
#         if not (z_bottom + margin < nozzle_z_abs < z_top - margin):
#             raise ValueError(
#                 f"nozzle_z_abs={nozzle_z_abs:.4f} m clips the diagrid face. "
#                 f"Must be in ({z_bottom + margin:.4f}, {z_top - margin:.4f})."
#             )

#         for theta in nozzle_boss_angles_deg:
#             rad = math.radians(theta)

#             # Point on the OUTER wall at this angle.
#             ox = radius * math.cos(rad)
#             oy = radius * math.sin(rad)

#             outward = cq.Vector( math.cos(rad),  math.sin(rad), 0.0)
#             inward  = cq.Vector(-math.cos(rad), -math.sin(rad), 0.0)
#             tangent = cq.Vector(-math.sin(rad),  math.cos(rad), 0.0)

#             # ── Base plane of the boss: INSET into the disc, NOT on the
#             #    outer wall. This is the key fix.
#             base_x = (radius - inset_safe) * math.cos(rad)
#             base_y = (radius - inset_safe) * math.sin(rad)
#             base_origin = cq.Vector(base_x, base_y, nozzle_z_abs)

#             boss_plane = cq.Plane(
#                 origin=base_origin, xDir=tangent, normal=outward
#             )
#             boss = (cq.Workplane(boss_plane)
#                     .circle(nozzle_r_boss)
#                     .extrude(boss_total_height))

#             # Union the boss with the disc — boolean union now produces a
#             # clean, gap-free junction because the boss base is embedded
#             # inside the cylindrical wall.
#             disc = disc.union(boss)

#             # ── Bore through the boss and into the diagrid body.
#             #    Start the bore at the OUTER face of the boss (i.e. at
#             #    `nozzle_boss_height` past the outer wall) and cut inward.
#             outer_face_origin = cq.Vector(
#                 (radius + nozzle_boss_height) * math.cos(rad),
#                 (radius + nozzle_boss_height) * math.sin(rad),
#                 nozzle_z_abs,
#             )
#             bore_plane = cq.Plane(
#                 origin=outer_face_origin, xDir=tangent, normal=inward
#             )
#             bore = (cq.Workplane(bore_plane)
#                     .circle(nozzle_r_bore)
#                     .extrude(nozzle_boss_height + nozzle_depth))
#             disc = disc.cut(bore)

#     return disc.clean()


# # ──────────────────────────────────────────────────────────────────────
# # Demo
# # ──────────────────────────────────────────────────────────────────────
# if __name__ == "__main__":
#     from ocp_vscode import show

#     # Reactor / pump geometry — kept in sync with centerline_assembly.py
#     _RPV_INNER_D = 8.91;  _RPV_WALL_T = 0.05
#     _TORI_Rc = 5.245;     _TORI_rk = 0.379
#     _od = _RPV_INNER_D + 2 * _RPV_WALL_T
#     _r  = _od / 2
#     _xk = _r - _TORI_rk
#     _zc = math.sqrt((_TORI_Rc - _TORI_rk)**2 - _xk**2)
#     _HEAD_BOTTOM_Z    = _zc - _TORI_Rc
#     _PUMP_Z_BOTTOM    = _HEAD_BOTTOM_Z + 2.562
#     _DIAGRID_Z_BOTTOM = -1.702 + 1.242
#     _pump_r           = 3.369

#     # Compute the azimuthal half-angle φ between the pump centerline and
#     # each nozzle mouth (same calculation as in centerline_assembly.py).
#     arc_rad   = math.radians(105.0)
#     L_leg     = 0.600; L_inlet = 0.050; R_bend = 0.460
#     barrel_r  = 1.350 / 2
#     overshoot = 0.040 * 1.5
#     ex  = R_bend * (1.0 - math.cos(arc_rad)) + L_leg * math.sin(arc_rad)
#     ey  = L_inlet + R_bend * math.sin(arc_rad) + L_leg * math.cos(arc_rad)
#     lx  = ey + (barrel_r - overshoot)
#     ly  = -ex
#     phi = math.degrees(math.atan2(lx, _pump_r + ly))

#     boss_angles = []
#     for a in [60.0, 180.0, 300.0]:
#         boss_angles.append(a - phi)
#         boss_angles.append(a + phi)

#     diagrid = create_diagrid(
#         diameter               = 4.660,
#         thickness              = 1.050,
#         z_bottom               = _DIAGRID_Z_BOTTOM,
#         nozzle_boss_angles_deg = boss_angles,
#         nozzle_z_abs           = _PUMP_Z_BOTTOM + 0.350,
#         nozzle_r_bore          = 0.230,
#         nozzle_depth           = 0.300,
#         nozzle_r_boss          = 0.301,
#         # Protrusion = (pump nozzle mouth radius) − (diagrid radius)
#         #            = 2.4075 − 2.3300 = 0.0775 m
#         nozzle_boss_height     = 0.0775,
#     )
#     show(diagrid)











































# """
# Parametric diagrid — ZLP V6.
# """

# from __future__ import annotations
# import math
# import cadquery as cq


# def create_diagrid(
#     diameter:               float,
#     thickness:              float,
#     z_bottom:               float = 0.0,
#     nozzle_boss_angles_deg: list[float] | None = None,
#     nozzle_z_abs:           float | None = None,
#     nozzle_r_bore:          float = 0.230,
#     nozzle_depth:           float = 0.300,
#     nozzle_r_boss:          float = 0.301,
#     nozzle_boss_height:     float = 0.080,
# ) -> cq.Workplane:

#     if diameter <= 0 or thickness <= 0:
#         raise ValueError("diameter and thickness must be > 0")

#     radius = diameter / 2.0
#     z_top  = z_bottom + thickness

#     disc = (cq.Workplane("XY")
#             .workplane(offset=z_bottom)
#             .circle(radius)
#             .extrude(thickness))

#     if nozzle_boss_angles_deg and nozzle_z_abs is not None:

#         margin = nozzle_r_boss
#         if not (z_bottom + margin < nozzle_z_abs < z_top - margin):
#             raise ValueError(
#                 f"nozzle_z_abs={nozzle_z_abs:.4f} m clips the diagrid face. "
#                 f"Must be in ({z_bottom + margin:.4f}, {z_top - margin:.4f})."
#             )

#         for theta in nozzle_boss_angles_deg:
#             rad = math.radians(theta)

#             ox = radius * math.cos(rad)
#             oy = radius * math.sin(rad)

#             outward = cq.Vector( math.cos(rad),  math.sin(rad), 0.0)
#             inward  = cq.Vector(-math.cos(rad), -math.sin(rad), 0.0)
#             tangent = cq.Vector(-math.sin(rad),  math.cos(rad), 0.0)
#             origin  = cq.Vector(ox, oy, nozzle_z_abs)

#             # Boss protruding outward
#             boss_plane = cq.Plane(origin=origin, xDir=tangent, normal=outward)
#             boss = (cq.Workplane(boss_plane)
#                     .circle(nozzle_r_boss)
#                     .extrude(nozzle_boss_height))
#             disc = disc.union(boss)

#             # Bore inward through boss and into diagrid body
#             bore_origin = origin + outward.multiply(nozzle_boss_height)
#             bore_plane  = cq.Plane(origin=bore_origin, xDir=tangent, normal=inward)
#             bore = (cq.Workplane(bore_plane)
#                     .circle(nozzle_r_bore)
#                     .extrude(nozzle_boss_height + nozzle_depth))
#             disc = disc.cut(bore)

#     return disc.clean()


# if __name__ == "__main__":
#     from ocp_vscode import show

#     _RPV_INNER_D = 8.91;  _RPV_WALL_T = 0.05
#     _TORI_Rc = 5.245;     _TORI_rk = 0.379
#     _od = _RPV_INNER_D + 2 * _RPV_WALL_T
#     _r  = _od / 2
#     _xk = _r - _TORI_rk
#     _zc = math.sqrt((_TORI_Rc - _TORI_rk)**2 - _xk**2)
#     _HEAD_BOTTOM_Z    = _zc - _TORI_Rc
#     _PUMP_Z_BOTTOM    = _HEAD_BOTTOM_Z + 2.562
#     _DIAGRID_Z_BOTTOM = -1.702 + 1.242
#     _pump_r           = 3.369

#     arc_rad   = math.radians(105.0)
#     L_leg     = 0.600; L_inlet = 0.050; R_bend = 0.460
#     barrel_r  = 1.350 / 2
#     overshoot = 0.040 * 1.5
#     ex  = R_bend * (1.0 - math.cos(arc_rad)) + L_leg * math.sin(arc_rad)
#     ey  = L_inlet + R_bend * math.sin(arc_rad) + L_leg * math.cos(arc_rad)
#     lx  = ey + (barrel_r - overshoot)
#     ly  = -ex
#     phi = math.degrees(math.atan2(lx, _pump_r + ly))

#     boss_angles = []
#     for a in [60.0, 180.0, 300.0]:
#         boss_angles.append(a - phi)
#         boss_angles.append(a + phi)

#     diagrid = create_diagrid(
#         diameter               = 4.660,
#         thickness              = 1.050,
#         z_bottom               = _DIAGRID_Z_BOTTOM,
#         nozzle_boss_angles_deg = boss_angles,
#         nozzle_z_abs           = _PUMP_Z_BOTTOM + 0.350,
#         nozzle_r_bore          = 0.230,
#         nozzle_depth           = 0.300,
#         nozzle_r_boss          = 0.301,
#         nozzle_boss_height     = 0.080,
#     )
#     show(diagrid)


















































# """
# Parametric diagrid — ZLP V7.
# Bore axes now follow the actual nozzle exit direction (not purely radial).
# """

# from __future__ import annotations
# import math
# import cadquery as cq


# def create_diagrid(
#     diameter:           float,
#     thickness:          float,
#     z_bottom:           float = 0.0,
#     nozzle_data:        list  | None = None,   # list of dicts: boss_angle_deg, bore_dx, bore_dy
#     nozzle_z_abs:       float | None = None,
#     nozzle_r_bore:      float = 0.230,
#     nozzle_depth:       float = 0.300,
#     nozzle_r_boss:      float = 0.301,
#     nozzle_boss_height: float = 0.080,
# ) -> cq.Workplane:

#     if diameter <= 0 or thickness <= 0:
#         raise ValueError("diameter and thickness must be > 0")

#     radius = diameter / 2.0
#     z_top  = z_bottom + thickness

#     disc = (cq.Workplane("XY")
#             .workplane(offset=z_bottom)
#             .circle(radius)
#             .extrude(thickness))

#     if nozzle_data and nozzle_z_abs is not None:

#         margin = nozzle_r_boss
#         if not (z_bottom + margin < nozzle_z_abs < z_top - margin):
#             raise ValueError(
#                 f"nozzle_z_abs={nozzle_z_abs:.4f} m clips the diagrid face. "
#                 f"Must be in ({z_bottom + margin:.4f}, {z_top - margin:.4f})."
#             )

#         for nozzle in nozzle_data:
#             theta   = math.radians(nozzle['boss_angle_deg'])
#             bore_dx = nozzle['bore_dx']   # inward unit vector x
#             bore_dy = nozzle['bore_dy']   # inward unit vector y

#             # Boss centre on diagrid outer curved surface
#             ox = radius * math.cos(theta)
#             oy = radius * math.sin(theta)

#             # Outward = opposite of bore (boss protrudes toward pump)
#             outward = cq.Vector(-bore_dx, -bore_dy, 0.0)
#             inward  = cq.Vector( bore_dx,  bore_dy, 0.0)
#             tangent = cq.Vector(-bore_dy,  bore_dx, 0.0)  # perpendicular to bore in XY

#             origin = cq.Vector(ox, oy, nozzle_z_abs)

#             # 1. Union cylindrical boss protruding outward toward pump
#             boss_plane = cq.Plane(origin=origin, xDir=tangent, normal=outward)
#             boss = (cq.Workplane(boss_plane)
#                     .circle(nozzle_r_boss)
#                     .extrude(nozzle_boss_height))
#             disc = disc.union(boss)

#             # 2. Cut bore from boss outer face inward, aligned with nozzle axis
#             bore_origin = origin + outward.multiply(nozzle_boss_height)
#             bore_plane  = cq.Plane(origin=bore_origin, xDir=tangent, normal=inward)
#             bore = (cq.Workplane(bore_plane)
#                     .circle(nozzle_r_bore)
#                     .extrude(nozzle_boss_height + nozzle_depth))
#             disc = disc.cut(bore)

#     return disc.clean()


# if __name__ == "__main__":
#     import math
#     from ocp_vscode import show

#     _RPV_INNER_D = 8.91;  _RPV_WALL_T = 0.05
#     _TORI_Rc = 5.245;     _TORI_rk = 0.379
#     _od = _RPV_INNER_D + 2 * _RPV_WALL_T;  _r = _od / 2
#     _xk = _r - _TORI_rk
#     _zc = math.sqrt((_TORI_Rc - _TORI_rk)**2 - _xk**2)
#     _HEAD_BOTTOM_Z    = _zc - _TORI_Rc
#     _PUMP_Z_BOTTOM    = _HEAD_BOTTOM_Z + 2.562
#     _DIAGRID_Z_BOTTOM = -1.702 + 1.242
#     _pump_r           = 3.369

#     arc_rad = math.radians(105.0)
#     barrel_r = 1.350/2; overshoot = 0.040*1.5
#     ex = 0.460*(1-math.cos(arc_rad)) + 0.600*math.sin(arc_rad)
#     ey = 0.050 + 0.460*math.sin(arc_rad) + 0.600*math.cos(arc_rad)
#     lx = ey + (barrel_r - overshoot)
#     ly = -ex
#     exit_r = ( math.cos(arc_rad), -math.sin(arc_rad))
#     exit_l = (-math.cos(arc_rad), -math.sin(arc_rad))

#     nozzle_data = []
#     for pump_angle_deg in [60.0, 180.0, 300.0]:
#         a = math.radians(pump_angle_deg)
#         sa, ca = math.sin(a), math.cos(a)
#         for (tlx, tly), (edx, edy) in [(( lx, ly), exit_r), ((-lx, ly), exit_l)]:
#             rx = tlx*sa + tly*ca;  ry = -tlx*ca + tly*sa
#             tip_x = _pump_r*ca + rx;  tip_y = _pump_r*sa + ry
#             nozzle_data.append({
#                 'boss_angle_deg': math.degrees(math.atan2(tip_y, tip_x)),
#                 'bore_dx': edx*sa + edy*ca,
#                 'bore_dy': -edx*ca + edy*sa,
#             })

#     diagrid = create_diagrid(
#         diameter           = 4.660,
#         thickness          = 1.050,
#         z_bottom           = _DIAGRID_Z_BOTTOM,
#         nozzle_data        = nozzle_data,
#         nozzle_z_abs       = _PUMP_Z_BOTTOM + 0.350,
#         nozzle_r_bore      = 0.230,
#         nozzle_depth       = 0.300,
#         nozzle_r_boss      = 0.301,
#         nozzle_boss_height = 0.080,
#     )
#     show(diagrid)



















































# """
# Parametric diagrid — ZLP V3.

# The diagrid is a horizontal disc that supports the reactor core assemblies
# in an SFR. Geometrically it's a flat cylindrical plate.

# All dimensions must be provided explicitly.

# Single public function:
#     create_diagrid()  — returns a single cq.Workplane solid
# """

# from __future__ import annotations
# import cadquery as cq


# def create_diagrid(
#     diameter:  float,
#     thickness: float,
#     z_bottom:  float = 0.0,
# ) -> cq.Workplane:
#     """
#     Build the diagrid as a flat cylindrical disc.

#     Args:
#         diameter:  outer diameter [m]
#         thickness: axial thickness [m]
#         z_bottom:  z coordinate of the bottom face [m] (default 0)

#     Returns:
#         cq.Workplane — disc with bottom face at z_bottom
#     """
#     if diameter <= 0:
#         raise ValueError(f"diameter must be > 0, got {diameter}")
#     if thickness <= 0:
#         raise ValueError(f"thickness must be > 0, got {thickness}")

#     return (cq.Workplane("XY")
#             .workplane(offset=z_bottom)
#             .circle(diameter / 2.0)
#             .extrude(thickness))


# # ---------------------------------------------------------------------------
# # Example
# # ---------------------------------------------------------------------------

# if __name__ == "__main__":
#     from ocp_vscode import show
#     diagrid = create_diagrid(
#         diameter  = 4.660,
#         thickness = 1.050,
#         z_bottom  = 0.0,
#     )
#     show(diagrid)


