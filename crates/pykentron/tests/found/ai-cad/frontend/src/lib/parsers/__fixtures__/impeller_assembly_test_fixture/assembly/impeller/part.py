"""Impeller — centrifugal compressor wheel (6 backswept blades).

Coordinate convention (per constraints.md, authoritative):
  bore axis = Z through origin; Z = 0 at the backplate BACKFACE
  (the running-clearance datum); +Z toward the hub nose / inducer eye.
  The assembly resolver places this local Z=0 at global +0.8.

Build chain (PROVEN by the pre-launch geometry probe — do not experiment):
  1. Single revolve: spigot + standoff boss + backplate + hub.
  2. ONE backswept blade: thin vertical profile swept along a plan-view
     spline path, then intersected with a revolved blade envelope
     (hub-to-shroud meridional region) to shape leading/trailing edges.
  3. Circular-pattern the blade 6x about Z, union all onto the hub.
  4. Cut the O10 through-bore LAST (avoids union-across-hollow failures).

NO fillets, NO chamfers — deferred project-wide.
"""

import math
import cadquery as cq

# === PARAMETERS (locked interface dims from constraints.md / design_plan.md) ===
D_MB_BORE = 17.0            # spigot OD: seats in 6203 bearing bore (press-fit intent, nominal)
Z_SPIGOT_END = -14.8        # spigot end face (shaft-shoulder clamp face)
Z_SPIGOT_TOP = -0.8         # spigot/boss junction (14.0 long = 12 bearing + 2 protrusion)
D_STANDOFF_BOSS = 22.0      # bears on 6203 inner-race face ONLY (race land ~O22.5, max O23)
T_STANDOFF_BOSS = 0.8       # boss thickness IS the running clearance CLR_RUN = 0.8
D_BACKPLATE = 72.0          # exducer tip diameter
T_BACKPLATE = 5.0           # backplate thickness: Z = 0 ... 5
D_HUB_BASE = 20.0           # hub diameter at backplate front face (Z = 5)
D_HUB_NOSE = 18.0           # hub nose diameter (lock-nut seat face; do NOT go below 18)
Z_NOSE = 34.0               # hub nose face (lock nut lands here; grip length 48.8 total)
D_SHAFT_IMPELLER_SEAT = 10.0  # through-bore, slide fit, modeled nominal
D_EYE = 34.0                # inducer eye diameter (blade LE tip circle at nose end)
N_BLADES = 6
T_BLADE = 2.0               # constant blade thickness (>= 1.5 min wall)
H_BLADE_EXDUCER = 6.0       # blade height above backplate front face at O72 rim
BACKSWEEP_DEG = 35.0        # backsweep angle at exducer (window 30-40 deg)

# === DERIVED VALUES ===
R_SPIGOT = D_MB_BORE / 2.0          # 8.5
R_BOSS = D_STANDOFF_BOSS / 2.0      # 11.0
R_TIP = D_BACKPLATE / 2.0           # 36.0
R_HUB_BASE = D_HUB_BASE / 2.0       # 10.0
R_HUB_NOSE = D_HUB_NOSE / 2.0       # 9.0
R_BORE = D_SHAFT_IMPELLER_SEAT / 2.0  # 5.0
R_EYE = D_EYE / 2.0                 # 17.0
Z_BP_FRONT = T_BACKPLATE            # 5.0, backplate front face
Z_RIM_TOP = Z_BP_FRONT + H_BLADE_EXDUCER  # 11.0, blade tip at exducer

# Blade leading edge stops 1.5 below the nose face so the M10 lock nut's
# across-corners overhang (~r9.25 > hub r9.0) can never touch a blade.
LE_NOSE_SETBACK = 1.5
Z_BLADE_TOP = Z_NOSE - LE_NOSE_SETBACK  # 32.5

# Blade camber (plan view): theta(r) = K_SWEEP * (r - R_BLADE_ROOT)^2.
# Backsweep angle beta at radius r satisfies tan(beta) = r * dtheta/dr,
# so K_SWEEP is solved for BACKSWEEP_DEG at the exducer tip (r = R_TIP).
R_BLADE_ROOT = 8.5          # camber start: buried inside hub (hub r >= 9.0), clear of bore r5
R_PATH_END = R_TIP + 1.0    # overshoot past rim: envelope trims trailing edge cleanly
K_SWEEP = math.tan(math.radians(BACKSWEEP_DEG)) / (
    2.0 * R_TIP * (R_TIP - R_BLADE_ROOT)
)

# Oversized blade sheet (trimmed by the envelope): generous Z overshoot both ends
Z_SHEET_BOT = 4.0           # below backplate front (buried -> robust union)
Z_SHEET_TOP = 36.0          # above blade top trim plane
R_ENV_INNER = 7.0           # envelope inner radius: inside hub (9..10), outside bore (5)
Z_ENV_BOT = 4.5             # envelope floor: inside backplate thickness
# Defensive offset (fix cycle 1): the buried blade foot must NOT share the
# r=36 cylinder face with the backplate rim (coincident faces made OCC fuse
# emit a zero-volume shard -> 2 solids, invalid). The floor stops at r=35 and
# a 1x1 diagonal reaches the full rim radius only ABOVE the backplate front.
R_ENV_FLOOR_OUT = 35.0      # buried-floor outer radius (1.0 inside the rim)
Z_ENV_RIM = 5.5             # z where the envelope first reaches r=36


# === FEATURE FUNCTIONS ===

def revolved_hub_body():
    """Spigot + standoff boss + backplate + tapered hub as ONE revolve.

    Half-profile in (r, z), drawn on the XZ plane and revolved about Z.
    Hub contour: straight monotonic taper O20 -> O18 (constraints allow a
    simple line profile; the sketch's curved silhouette is the blade shroud,
    reproduced by blade_envelope(), not the hub itself).
    """
    profile_pts = [
        (0.0, Z_SPIGOT_END),          # axis, spigot end face
        (R_SPIGOT, Z_SPIGOT_END),     # spigot end face rim
        (R_SPIGOT, Z_SPIGOT_TOP),     # spigot OD up to boss underside
        (R_BOSS, Z_SPIGOT_TOP),       # boss underside (inner-race contact face)
        (R_BOSS, 0.0),                # boss OD, 0.8 thick
        (R_TIP, 0.0),                 # backplate BACKFACE = Z datum
        (R_TIP, Z_BP_FRONT),          # backplate rim, 5.0 thick
        (R_HUB_BASE, Z_BP_FRONT),     # backplate front, in to hub base O20
        (R_HUB_NOSE, Z_NOSE),         # hub taper to nose O18
        (0.0, Z_NOSE),                # nose face in to axis
    ]
    return (
        cq.Workplane("XZ")
        .polyline(profile_pts)
        .close()
        .revolve(360, (0, 0, 0), (0, 1, 0))
    )


def blade_envelope():
    """Revolved meridional region between hub and shroud — intersection tool.

    Outer boundary: rim cylinder r36 up to z11 (6.0 exducer height), then the
    shroud spline curving up to the eye radius r17 near the nose (this is the
    curved blade-tip silhouette from the user sketch). Inner boundary r7 is
    buried inside the hub. Floor z4.5 is buried inside the backplate.
    """
    shroud_pts = [
        (30.0, 13.0),
        (24.0, 17.5),
        (20.0, 24.0),
        (18.0, 29.0),
        (R_EYE, Z_BLADE_TOP),
    ]
    return (
        cq.Workplane("XZ")
        .moveTo(R_ENV_INNER, Z_ENV_BOT)
        .lineTo(R_ENV_FLOOR_OUT, Z_ENV_BOT)
        .lineTo(R_TIP, Z_ENV_RIM)
        .lineTo(R_TIP, Z_RIM_TOP)
        .spline(shroud_pts, includeCurrent=True)
        .lineTo(R_ENV_INNER, Z_BLADE_TOP)
        .close()
        .revolve(360, (0, 0, 0), (0, 1, 0))
    )


def backswept_blade(envelope):
    """ONE backswept blade: sweep a thin profile along a spline path (PROVEN).

    Path: plan-view camber spline theta = K_SWEEP*(r - root)^2 at Z=0
    (start tangent is purely radial, so the profile plane normal +X matches).
    Profile: T_BLADE x 32 vertical rectangle on the YZ plane at the path start.
    The swept sheet is translated up, then trimmed by the envelope intersect.
    """
    radii = [8.5, 14.0, 20.0, 26.0, 31.0, 34.0, R_PATH_END]
    path_pts = []
    for r in radii:
        th = K_SWEEP * (r - R_BLADE_ROOT) ** 2
        path_pts.append((r * math.cos(th), r * math.sin(th)))

    path = cq.Workplane("XY").spline(path_pts)

    sheet_h = Z_SHEET_TOP - Z_SHEET_BOT  # 32.0
    sheet = (
        cq.Workplane("YZ")
        .workplane(offset=R_BLADE_ROOT)   # profile plane at path start (x=8.5)
        .rect(T_BLADE, sheet_h)
        .sweep(path, isFrenet=True)
    )
    # Profile was centered at z=0: lift sheet so it spans Z_SHEET_BOT..Z_SHEET_TOP
    sheet = sheet.translate((0, 0, (Z_SHEET_TOP + Z_SHEET_BOT) / 2.0))

    return sheet.intersect(envelope)


def blade_array(body, blade):
    """Circular-pattern the blade N_BLADES x about Z and union onto the hub."""
    for i in range(N_BLADES):
        angle = i * (360.0 / N_BLADES)
        body = body.union(blade.rotate((0, 0, 0), (0, 0, 1), angle))
    return body


def shaft_through_bore(body):
    """O10 through-bore, cut LAST (proven order: union blades first, then bore)."""
    cutter = (
        cq.Workplane("XY")
        .workplane(offset=Z_SPIGOT_END - 1.0)
        .circle(R_BORE)
        .extrude((Z_NOSE - Z_SPIGOT_END) + 2.0)
    )
    return body.cut(cutter)


# === COMPOSE (reads like the manufacturing sequence) ===
body = revolved_hub_body()
envelope = blade_envelope()
blade = backswept_blade(envelope)
body = blade_array(body, blade)
body = shaft_through_bore(body)
result = body

# === NUMERIC SELF-CHECKS (printed for the executor log; vision is deferred) ===
_solids = result.solids().vals()
_shape = result.val()
_bb = _shape.BoundingBox()
print(f"CHECK solids={len(_solids)}")
print(f"CHECK valid={_shape.isValid()}")
print(
    "CHECK bbox="
    f"x[{_bb.xmin:.3f},{_bb.xmax:.3f}] "
    f"y[{_bb.ymin:.3f},{_bb.ymax:.3f}] "
    f"z[{_bb.zmin:.3f},{_bb.zmax:.3f}] "
    f"lens=({_bb.xlen:.3f},{_bb.ylen:.3f},{_bb.zlen:.3f})"
)
print(f"CHECK volume={_shape.Volume():.1f} mm^3")
