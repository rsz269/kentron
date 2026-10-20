"""Side rail for dual MacBook Pro 13-inch shelf stand (two-part with diagonal tenon).

Print 2 copies of each rail part (4 total: 2 front, 2 back).
Print 8 bridge bars (4 lower + 4 upper).
Assembly:
  1. Join front and back rail halves at 15-degree angle (diagonal tenon).
  2. Connect left and right rails with 8 bridges (slots on ledge surfaces).
  3. Place assembled frame in shelf (free-standing, not against walls).

Cross-section (front view, one rail shown):

    A----B
    |    |
    |    C-------D     upper ledge (upper MacBook sits here)
    |    |       |
    |    |       E
    |    |     /F      gusset (structural + print-friendly overhang)
    |    |   /
    |    G
    |    |
    |    H-------I     lower ledge / base top (lower MacBook sits here)
    |    |       |
    K-----------J      base bottom (vent slots cut here)

Frame assembly (top view):

    left rail                          right rail
    body|ledge|                        |ledge|body
    ┌───┤═════╪════════════════════════╪═════├───┐
    │   │     │ bridge ①  (Z=5)       │     │   │
    │   │     │ bridge ②  (Z=80)      │     │   │
    │   │     │    [joint zone]        │     │   │
    │   │     │ bridge ③  (Z=178)     │     │   │
    │   │     │ bridge ④  (Z=253)     │     │   │
    └───┤═════╪════════════════════════╪═════├───┘
         58mm  ←──── 130mm gap ────→   58mm
              ←────── 230mm inner ──────→

Recommended print: body wall flat on bed, 58mm print height, no supports needed.
"""

import math
from pathlib import Path

import cadquery as cq

# === MacBook Pro 13" ===
mb_width = 313.0  # mm - body width
mb_depth = 220.0  # mm - body depth
mb_thickness = 18.0  # mm - closed-lid thickness

# === Shelf interior (3 walls + ceiling, front open) ===
shelf_width = 385.0  # mm - internal width
shelf_depth = 261.0  # mm - internal depth
shelf_height = 135.0  # mm - internal height

# === Print constraints (Bambu Lab A1 Mini) ===
max_print = 180.0  # mm - build volume per axis

# === Design parameters ===
body_t = 8.0  # mm - body wall thickness (structural spine)
ledge_w = 50.0  # mm - ledge extending inward from body (MacBook support)
ledge_t = 8.0  # mm - ledge platform thickness (was 5, increased for joint strength)
foot_h = 10.0  # mm - base foot height (airflow channel under lower MacBook)
air_gap = 40.0  # mm - gap between MacBooks (ventilation + finger access for removal)
gusset_reach = 15.0  # mm - gusset extent along upper ledge underside
gusset_drop = 20.0  # mm - gusset vertical drop below upper ledge
rail_depth = 258.0  # mm - total rail depth (178 original + 80mm extension)
top_ext = 8.0  # mm - body wall extending above upper ledge top surface
fit_tol = 0.2  # mm - general sliding fit tolerance

# === Frame / bridge parameters ===
inner_gap = 230.0  # mm - distance between body wall inner faces (MacBook 220 + 5mm x2)
bridge_bar_w = 16.0  # mm - bridge bar width along Z (wider than tab for Z shoulder)
bridge_bar_h = 2.0  # mm - bridge bar height (Y, recessed flush into ledge surface)
bridge_tab_len = 15.0  # mm - tab length in X (how far tab extends onto ledge)
bridge_tab_w = 10.0  # mm - tab width in Z (narrower than bar, 3mm shoulder each side)
pocket_wall = 3.0  # mm - min wall thickness between pocket and ledge tip (X direction)
pocket_floor = 2.0  # mm - solid floor at pocket bottom (prevents falling through)

# === Cantilever snap-fit parameters ===
cant_t = 1.0  # mm - cantilever arm thickness (X direction, flex)
cant_clearance = 0.5  # mm - FDM assembly clearance per side
hook_overhang = 0.4  # mm - hook protrusion (X, deflection ~3.75% strain for PLA)
hook_h = 1.0  # mm - hook height (Y direction)
# End bridges: center = bar_w/2 so edges align with rail edges (Z=0, Z=258)
bridge_z_positions = [8.0, 80.0, 178.0, 250.0]  # Z centers

# === Diagonal tenon joint parameters ===
tenon_angle = 15.0  # degrees - tenon tilt from Z axis toward +X (locks X and Z axes)
tenon_length = 12.0  # mm - tenon length along its tilted axis
tenon_width = 40.0  # mm - tenon width in X direction
tenon_taper = 0.3  # mm - groove narrows by this much from entry to deep end (wedge)
tenon_tol_entry = 0.20  # mm - clearance per side at entry (loose, easy start)
tenon_tol_deep = 0.10  # mm - clearance per side at deep end (tight, friction lock)
tenon_h_base = 10.0  # mm - base tenon height (fits in 18mm base section)
tenon_h_upper = 6.0  # mm - upper tenon height (fits in 8mm upper ledge, was 4mm)

# === Ventilation slots (cut through base foot area) ===
vent_w = 8.0  # mm - slot width (67% open ratio)
vent_pitch = 12.0  # mm - slot center-to-center spacing
vent_margin_x = 8.0  # mm - inset from body inner face and ledge edge
vent_margin_z = 6.0  # mm - inset from front and back of rail
vent_bridge_z = 30.0  # mm - solid bridge width at joint (centered on z_split)

# === Body wall ventilation holes (hex grid pattern) ===
hole_d = 6.0  # mm - hole diameter
hole_pitch = 10.0  # mm - center-to-center hex grid spacing
hole_margin_y = 3.0  # mm - clearance from structural Y features
hole_margin_z = 6.0  # mm - clearance from Z edges
hole_tenon_clearance = 15.0  # mm - clearance around tenon joint zone

# === Airflow spacer ribs (on body wall inner face) ===
rib_depth = 2.5  # mm - rib protrusion from body wall inner face (+X direction)
rib_width = 3.0  # mm - rib thickness along Z direction
rib_count = 3  # number of ribs evenly spaced along Z

# === Derived dimensions ===
base_h = foot_h + ledge_t  # 18 mm - total base section height
upper_z = base_h + mb_thickness + air_gap  # 76 mm - bottom of upper ledge
rail_h = upper_z + ledge_t + top_ext  # 92 mm - total rail height
z_split = rail_depth / 2  # 129 mm - split point (equal halves)

# Bridge derived dimensions
bridge_span = inner_gap - 2 * ledge_w  # 130 mm - gap between ledge tips
bridge_total_len = bridge_span + 2 * bridge_tab_len  # 160 mm - total bar length

# Tenon extent in global coordinates
tenon_z_extent = tenon_length * math.cos(math.radians(tenon_angle))
tenon_x_extent = tenon_length * math.sin(math.radians(tenon_angle))

# === Validation ===
part_max_z = z_split + tenon_z_extent + (tenon_width / 2) * math.sin(math.radians(tenon_angle))
assert part_max_z <= max_print, f"Part depth {part_max_z:.1f} exceeds build volume"
assert rail_h <= max_print, f"Rail height {rail_h} exceeds build volume {max_print}"
assert body_t + ledge_w <= max_print, "Rail width exceeds build volume"
assert rail_h < shelf_height, f"Rail height {rail_h} exceeds shelf height"
assert upper_z - gusset_drop > base_h, "Gusset overlaps base section"
assert bridge_total_len <= max_print, f"Bridge {bridge_total_len:.0f} exceeds build volume"

# === Build cross-section profile ===
# Coordinate system: X = width (outer wall -> ledge), Y = height (base -> top)
# Extrude along Z = depth (back -> front)

pts = [
    (0, rail_h),  # A: outer top
    (body_t, rail_h),  # B: inner top
    (body_t, upper_z + ledge_t),  # C: body at upper ledge top
    (body_t + ledge_w, upper_z + ledge_t),  # D: upper ledge tip, top
    (body_t + ledge_w, upper_z),  # E: upper ledge tip, bottom
    (body_t + gusset_reach, upper_z),  # F: gusset start on ledge underside
    (body_t, upper_z - gusset_drop),  # G: gusset end on body
    (body_t, base_h),  # H: body at base top / lower ledge surface
    (body_t + ledge_w, base_h),  # I: base top, ledge tip
    (body_t + ledge_w, 0),  # J: base bottom, outer edge
    (0, 0),  # K: outer bottom
]

full_rail = cq.Workplane("XY").moveTo(*pts[0]).polyline(pts[1:]).close().extrude(rail_depth)

# === Cut ventilation slots through foot area of base ===
# Slots split into front/back segments with solid bridge at joint for strength.
slot_x_start = body_t + vent_margin_x
slot_x_end = body_t + ledge_w - vent_margin_x

slot_front_len = z_split - vent_bridge_z / 2 - vent_margin_z
slot_front_ctr = vent_margin_z + slot_front_len / 2
slot_back_len = rail_depth - vent_margin_z - (z_split + vent_bridge_z / 2)
slot_back_ctr = z_split + vent_bridge_z / 2 + slot_back_len / 2

cx = slot_x_start
while cx + vent_w / 2 <= slot_x_end:
    for s_ctr, s_len in [
        (slot_front_ctr, slot_front_len),
        (slot_back_ctr, slot_back_len),
    ]:
        full_rail = full_rail.cut(
            cq.Workplane("XY")
            .transformed(offset=(cx, foot_h / 2, s_ctr))
            .box(vent_w, foot_h + 0.2, s_len)
        )
    cx += vent_pitch

# === Cut bridge snap-fit pockets in ledges ===
# Blind pocket with undercut groove for cantilever snap-fit.
# Upper portion: narrow (arm passage with deflection clearance).
# Lower portion: wider on inner side (hook engagement groove).
# Solid floor at bottom prevents falling through.
# Assembly: push bridge down from above, arm deflects, hook snaps into groove.

pocket_depth = ledge_t - pocket_floor  # 6mm deep (2mm floor remains)
pocket_z = bridge_tab_w + 2 * cant_clearance  # pocket width in Z

# Pocket X dimensions (at ledge tip area)
pocket_x_main = cant_t + 2 * cant_clearance  # narrow upper: arm + clearance
pocket_x_groove = hook_overhang + cant_clearance  # groove: extra width for hook
pocket_x_total = pocket_x_main + pocket_x_groove  # total at hook level
pocket_x_outer = body_t + ledge_w - pocket_wall  # outer edge position

ledge_y_tops = [
    base_h,  # lower ledge surface: Y=18
    upper_z + ledge_t,  # upper ledge surface: Y=84
]

bar_recess_z = bridge_bar_w + 2 * cant_clearance  # recess width in Z for bar

for y_top in ledge_y_tops:
    for bz in bridge_z_positions:
        y_floor = y_top - pocket_depth
        # Bar recess: shallow cut from ledge tip inward, so bar sits flush.
        # Must extend to the tip face (X=body_t+ledge_w) for bar entry.
        recess_x_center = body_t + ledge_w - bridge_tab_len / 2
        full_rail = full_rail.cut(
            cq.Workplane("XY")
            .transformed(
                offset=(
                    recess_x_center,
                    y_top - bridge_bar_h / 2,
                    bz,
                )
            )
            .box(bridge_tab_len + 0.1, bridge_bar_h + 0.1, bar_recess_z)
        )
        # Main pocket (narrow, full depth): arm passage
        full_rail = full_rail.cut(
            cq.Workplane("XY")
            .transformed(
                offset=(
                    pocket_x_outer - pocket_x_main / 2,
                    (y_top + y_floor) / 2,
                    bz,
                )
            )
            .box(pocket_x_main, pocket_depth + 0.1, pocket_z)
        )
        # Hook groove on OUTER wall (tip side) for outward-force resistance
        groove_height = hook_h + cant_clearance
        full_rail = full_rail.cut(
            cq.Workplane("XY")
            .transformed(
                offset=(
                    pocket_x_outer + pocket_x_groove / 2,
                    y_floor + groove_height / 2,
                    bz,
                )
            )
            .box(pocket_x_groove, groove_height, pocket_z)
        )

# === Cut body wall ventilation holes (hex grid) ===
hole_r = hole_d / 2
hex_row_spacing = hole_pitch * math.sqrt(3) / 2

y_exclusions = [
    (0, hole_margin_y),
    (base_h - hole_margin_y, base_h + hole_margin_y),
    (
        upper_z - gusset_drop - hole_margin_y,
        upper_z + ledge_t + hole_margin_y,
    ),
    (rail_h - hole_margin_y, rail_h),
]

z_hole_min = hole_margin_z + hole_r
z_hole_max = rail_depth - hole_margin_z - hole_r
z_excl_min = z_split - hole_tenon_clearance
z_excl_max = z_split + hole_tenon_clearance

hole_positions = []
y = hole_margin_y + hole_r
row = 0
while y <= rail_h - hole_margin_y - hole_r:
    if not any(ylo <= y <= yhi for ylo, yhi in y_exclusions):
        z_offset = (hole_pitch / 2) if (row % 2 == 1) else 0
        z = z_hole_min + z_offset
        while z <= z_hole_max:
            if not (z_excl_min <= z <= z_excl_max):
                hole_positions.append((y, z))
            z += hole_pitch
    y += hex_row_spacing
    row += 1

if hole_positions:
    body_holes = (
        cq.Workplane("YZ")
        .transformed(offset=(0, 0, -0.1))
        .pushPoints(hole_positions)
        .circle(hole_r)
        .extrude(body_t + 0.2)
    )
    full_rail = full_rail.cut(body_holes)

# === Airflow spacer ribs on body wall inner face ===
rib_y_start = base_h
rib_y_end = upper_z + ledge_t
rib_height = rib_y_end - rib_y_start
rib_z_spacing = rail_depth / (rib_count + 1)

for i in range(rib_count):
    rib_z = rib_z_spacing * (i + 1)
    full_rail = full_rail.union(
        cq.Workplane("XY")
        .transformed(
            offset=(
                body_t + rib_depth / 2,
                rib_y_start + rib_height / 2,
                rib_z,
            )
        )
        .box(rib_depth, rib_height, rib_width)
    )

# === Diagonal tenon joint ===
margin = 2.0
bb_w = body_t + ledge_w + 2 * margin
bb_h = rail_h + 2 * margin
bb_cx = (body_t + ledge_w) / 2
bb_cy = rail_h / 2

tenon_cx = (body_t + ledge_w) / 2

tenon_positions = [
    (tenon_cy_base := base_h / 2, tenon_h_base),
    (tenon_cy_upper := upper_z + ledge_t / 2, tenon_h_upper),
]


def make_tenon_pair(cy, th):
    """Create tongue (undersized) and groove (nominal) for one tenon position."""
    tongue = (
        cq.Workplane("XY")
        .transformed(offset=(tenon_cx, cy, z_split), rotate=(0, tenon_angle, 0))
        .rect(tenon_width - 2 * tenon_tol_entry, th - 2 * tenon_tol_entry)
        .workplane(offset=tenon_length - tenon_tol_entry)
        .rect(
            tenon_width - tenon_taper - 2 * tenon_tol_deep,
            th - tenon_taper - 2 * tenon_tol_deep,
        )
        .loft()
    )
    groove = (
        cq.Workplane("XY")
        .transformed(offset=(tenon_cx, cy, z_split), rotate=(0, tenon_angle, 0))
        .rect(tenon_width, th)
        .workplane(offset=tenon_length)
        .rect(tenon_width - tenon_taper, th - tenon_taper)
        .loft()
    )
    return tongue, groove


# === Split rail into two parts ===
front_box = (
    cq.Workplane("XY").transformed(offset=(bb_cx, bb_cy, z_split / 2)).box(bb_w, bb_h, z_split)
)
back_box = (
    cq.Workplane("XY")
    .transformed(offset=(bb_cx, bb_cy, z_split + (rail_depth - z_split) / 2))
    .box(bb_w, bb_h, rail_depth - z_split)
)

front_cutter = front_box
back_cutter = back_box

for cy, th in tenon_positions:
    tongue, groove = make_tenon_pair(cy, th)
    front_cutter = front_cutter.union(tongue)
    back_cutter = back_cutter.cut(groove)

part_front = full_rail.intersect(front_cutter)
part_back = full_rail.intersect(back_cutter)

# === Build bridge bar ===
# Bar with L-shaped hooks at each end that wrap around the ledge tips.
# The hooks provide mechanical locking against X-axis forces.
#
# Cross-section at each end (X-Y, looking along Z):
#
#     ═══════════ bar (on ledge surface)
#          │
#          │  vertical (down ledge outer face)
#          │
#          └── hook (under ledge, X-lock)
#
# Bridge X coords: 0 = left hook outer edge, fl = right hook outer edge
# Bridge Y coords: 0 = ledge surface level, positive = up

# Bridge: bar with cantilever snap-fit arms at each end.
# Built as a single XZ profile extruded along Y, then cut to create arms.
# Arms narrower than bar (Z direction) for shoulder-based Z lock.
# Assembly: push bridge down onto ledge, arms deflect and snap.

fl = bridge_total_len
arm_z = bridge_tab_w - 2 * cant_clearance  # arm Z width (narrower than bar)
arm_depth = pocket_depth - cant_clearance  # arm Y height

# Arm/hook positions
arm_inset = pocket_wall + cant_clearance  # 3.5mm from bridge edge
left_arm_cx = arm_inset + cant_t / 2  # 4.0
right_arm_cx = fl - arm_inset - cant_t / 2  # fl - 4.0

# Build bridge as XZ profile extruded along Y (arm_z width).
# Profile includes: bar + left arm + left hook + right arm + right hook.
# Then widen the bar portion to bridge_bar_w using a separate union.

# Step 1: Central bar as a simple box (full width in Y)
bridge = (
    cq.Workplane("XY")
    .transformed(offset=(fl / 2, 0, bridge_bar_h / 2))
    .box(fl, bridge_bar_w, bridge_bar_h)
)

# Step 2: Build arms+hooks as a single XZ profile, extruded narrow (arm_z)
# This avoids union issues by creating one solid per arm.
for i in range(2):
    direction = -1 if i == 0 else 1  # outward direction
    arm_cx = left_arm_cx if i == 0 else right_arm_cx
    outer_x = arm_cx + direction * (cant_t / 2)

    # Profile: arm rectangle + hook triangle, as one connected polygon
    # Arm extends from Z=+1 (overlap into bar) to Z=-arm_depth
    # Hook wedge at bottom extends outward
    arm_pts = (
        [
            (arm_cx - cant_t / 2, 1.0),  # arm top-inner (into bar)
            (arm_cx + cant_t / 2, 1.0),  # arm top-outer (into bar)
            (arm_cx + cant_t / 2, -arm_depth + hook_h),  # arm bottom, before hook
            # Hook ramp (wedge "4" shape)
            (outer_x + direction * hook_overhang, -arm_depth),  # hook tip
            (arm_cx - cant_t / 2, -arm_depth),  # hook base / arm inner bottom
        ]
        if direction == -1
        else [
            (arm_cx + cant_t / 2, 1.0),  # arm top-outer
            (arm_cx - cant_t / 2, 1.0),  # arm top-inner (into bar)
            (arm_cx - cant_t / 2, -arm_depth + hook_h),  # arm bottom, before hook
            # Hook ramp (wedge "4" shape)
            (outer_x + direction * hook_overhang, -arm_depth),  # hook tip
            (arm_cx + cant_t / 2, -arm_depth),  # hook base / arm outer bottom
        ]
    )

    arm_solid = (
        cq.Workplane("XZ")
        .transformed(offset=(0, 0, -arm_z / 2))
        .moveTo(*arm_pts[0])
        .polyline(arm_pts[1:])
        .close()
        .extrude(arm_z)
    )
    bridge = bridge.union(arm_solid)

# === Export ===
output_dir = Path(__file__).parent / "output"
output_dir.mkdir(exist_ok=True)

front_path = output_dir / "macbook-stand-front.stl"
back_path = output_dir / "macbook-stand-back.stl"
bridge_path = output_dir / "macbook-stand-bridge.stl"

part_front.export(str(front_path))
part_back.export(str(back_path))
bridge.export(str(bridge_path))

try:
    show_object(part_front, name="front")  # type: ignore[name-defined]
    show_object(part_back, name="back")  # type: ignore[name-defined]
    show_object(bridge, name="bridge")  # type: ignore[name-defined]
except NameError:
    pass

print(f"Total rail: {body_t + ledge_w:.0f} x {rail_depth:.0f} x {rail_h:.0f} mm")
print(f"Split at Z={z_split:.0f} mm")
print(f"Base tenon: {tenon_width:.0f}x{tenon_h_base:.0f} mm at Y={tenon_cy_base:.1f}")
print(f"Upper tenon: {tenon_width:.0f}x{tenon_h_upper:.0f} mm at Y={tenon_cy_upper:.1f}")
print(f"Tenon L={tenon_length:.0f} mm, angle={tenon_angle:.0f} deg")
print(f"Frame inner gap: {inner_gap:.0f} mm (MacBook {mb_depth:.0f} + clearance)")
print(f"Bridge: {bridge_total_len:.0f} mm (span {bridge_span:.0f} + tabs {bridge_tab_len:.0f}x2)")
print(f"Bridge positions: {bridge_z_positions} ({len(bridge_z_positions)} per level, x2 levels)")
print(f"Body wall holes: {len(hole_positions)} x dia {hole_d:.0f} mm (hex grid)")
print(f"Airflow ribs: {rib_count} x {rib_depth:.1f}mm deep (air gap for exhaust)")
vent_open_pct = vent_w / vent_pitch * 100
print(f"Base vent slots: {vent_w:.0f}mm / {vent_pitch:.0f}mm pitch ({vent_open_pct:.0f}% open)")
print(f"Each part max depth: {part_max_z:.0f} mm (limit: {max_print:.0f} mm)")
print(f"Exported: {front_path}")
print(f"Exported: {back_path}")
print(f"Exported: {bridge_path} (print 8 copies)")
