import cadquery as cq

# =============================================================================
# lid — half-cylinder shell lid for the treasure chest (FDM printed)
#
# Local frame (per constraints.md): centered on XY, flat rim face at Z=0,
# dome in +Z, half-cylinder axis parallel to X at (Y=0, Z=0). Rear = -Y.
# In the assembly, local Z=0 == global Z=80 (base rim plane), so the local
# hinge axis (Y=-66, Z=0) == global (Y=-66, Z=80).
# NO fillets (deferred project-wide).
#
# PROPOSAL (assembly_resolver, CONFLICT-001): identical to part.py except for
# the added `lug_underhang_trim` feature. The lugs embed to Y=-57 (inner dome
# radius) for a solid union — correct for local Z>=0 where the 3 mm dome wall
# exists, but below the rim (local Z=-4..0) there is no dome material, and the
# bare Y=-60..-57 slab lands inside the base_box rear wall (Y=-60..-57,
# Z=0..80) when the lid is placed at global Z=+80: a 15 x 3 x 4 mm
# interpenetration per lug. The trim restores the constraints.md geometry
# below the rim (tab spans Y=-74..-60, flush on the base outer wall).
# =============================================================================

# === PARAMETERS (mm) ===
# Primary — from constraints.md
dome_outer_radius = 60.0     # half-cylinder outer radius
wall_thickness = 3.0         # shell wall AND end-cap thickness
lid_length = 180.0           # overall length along X (matches base footprint)
lug_width_x = 15.0           # hinge lug width along X
lug_center_x = 37.0          # lug centers at X = +37 / -37 (spans |X| 29.5..44.5)
lug_rear_y = -74.0           # lug rear face (14 mm rearward of dome at Y=-60)
lug_half_height_z = 4.0      # lug spans Z = -4..+4 (4 mm below the rim plane)
lug_bore_dia = 3.4           # clearance bore for ISO 2338 dowel pin d3 x 35
                             # (~0.4 mm running clearance — pin ROTATES in lug;
                             # nominal 3.4, NO FDM compensation per sourcing)
hinge_axis_y = -66.0         # bore axis, part-local
hinge_axis_z = 0.0           # bore axis on the rim plane (local Z=0)

# Derived
dome_inner_radius = dome_outer_radius - wall_thickness   # 57.0
half_length = lid_length / 2.0                           # 90.0
cavity_length = lid_length - 2.0 * wall_thickness        # 174.0 -> 3 mm end caps
# Lug body embeds to the INNER radius so the union with the curved shell is
# fully solid (a lug ending exactly at Y=-60 would only touch the dome outer
# surface along the Z=0 tangent line -> degenerate contact). Embedding to
# Y=-57 keeps the lug inside the 3 mm wall material without intruding into
# the interior cavity (the inner surface curves inward for Z>0). BELOW the
# rim plane (Z<0) there is no wall to embed into — that portion is trimmed
# back to Y=-60 by lug_underhang_trim (CONFLICT-001).
lug_front_y = -dome_inner_radius                         # -57.0
lug_span_y = lug_front_y - lug_rear_y                    # 17.0
lug_center_y = (lug_rear_y + lug_front_y) / 2.0          # -65.5
dome_rear_y = -dome_outer_radius                         # -60.0 — base wall outer face


# === FEATURE FUNCTIONS ===
def dome_blank():
    """Solid half-cylinder R60 x 180: flat face on Z=0, axis along X through origin."""
    profile = (
        cq.Workplane("YZ")  # local (u, v) = (global Y, global Z), normal = +X
        .moveTo(-dome_outer_radius, 0.0)
        .threePointArc((0.0, dome_outer_radius), (dome_outer_radius, 0.0))
        .close()
    )
    return profile.extrude(half_length, both=True)  # X = -90 .. +90


def interior_cavity(body):
    """Hollow the dome: R57 x 174 semicylindrical cavity.

    Leaves a 3 mm curved wall, 3 mm flat semicircular end caps at X=+/-87..90,
    and opens the flat rim face (local Z=0). The cutter carries a 1 mm skirt
    below Z=0 so the boolean at the open rim is clean (no coplanar-face cut).
    """
    skirt = 1.0
    cutter = (
        cq.Workplane("YZ")
        .moveTo(-dome_inner_radius, -skirt)
        .lineTo(dome_inner_radius, -skirt)
        .lineTo(dome_inner_radius, 0.0)
        .threePointArc((0.0, dome_inner_radius), (-dome_inner_radius, 0.0))
        .close()
        .extrude(cavity_length / 2.0, both=True)  # X = -87 .. +87
    )
    return body.cut(cutter)


def hinge_lugs(body):
    """Two integral rear hinge lugs: 15 wide (X), Y=-74..-57, Z=-4..+4.

    Visible projection is Y=-74..-60 and 4 mm below the rim, per constraints;
    the extra Y=-60..-57 portion is buried in the shell wall for a solid merge
    (valid only for Z>=0 — the Z<0 part of that slab is trimmed next).
    """
    lug_height_z = 2.0 * lug_half_height_z  # 8.0 -> >=2.3 mm plastic around the bore
    for cx in (+lug_center_x, -lug_center_x):
        lug = (
            cq.Workplane("XY")
            .box(lug_width_x, lug_span_y, lug_height_z)  # centered on origin
            .translate((cx, lug_center_y, 0.0))
        )
        body = body.union(lug)
    return body


def lug_underhang_trim(body):
    """Trim the lug embedment slab BELOW the rim plane back to Y=-60 (CONFLICT-001).

    The dome shell exists only for Z>=0, so the Y=-60..-57 embedment is bare
    below the rim and would interpenetrate the base_box rear wall (which spans
    Y=-60..-57 up to Z=80) when the lid is closed. One cut box spanning
    Y=-60..-50, Z=-8..0, full X width removes exactly that slab from both lugs
    (the dome has no material at Z<0, so nothing else is touched). Below the
    rim the lugs now span Y=-74..-60 — flush contact with the base outer wall,
    exactly as constraints.md specifies.
    """
    trim = (
        cq.Workplane("XY")
        .box(lid_length + 10.0, 10.0, 2.0 * lug_half_height_z)
        .translate((0.0, dome_rear_y + 5.0, -lug_half_height_z))
    )  # X = -95..+95, Y = -60..-50, Z = -8..0
    return body.cut(trim)


def hinge_pin_bores(body):
    """D3.4 through-bores along X at (Y=-66, Z=0) — one cutter through both lugs.

    The axis lies rearward of the dome (dome rear extent is Y=-60), so a single
    full-length cylinder only intersects the two lugs.
    """
    cutter = (
        cq.Workplane("YZ")
        .center(hinge_axis_y, hinge_axis_z)
        .circle(lug_bore_dia / 2.0)
        .extrude(half_length + 5.0, both=True)  # X = -95 .. +95, safely through
    )
    return body.cut(cutter)


# === COMPOSE (reads like the build sequence) ===
body = dome_blank()
body = interior_cavity(body)
body = hinge_lugs(body)
body = lug_underhang_trim(body)
body = hinge_pin_bores(body)
result = body
