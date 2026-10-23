# Generated with YoCAD — https://yotown.com · 2026-07-18 · run 019f779b-c061-7712-b2f0-cc7384db3318
"""YoCAD 4D concept — "ARES-01" Mars cargo rover.

A compact pressurized central body (grounded root) carries an open, actively
leveled rear cargo deck, four large driven wheels on rotating wheel modules,
and four articulated robotic legs (deploy → hip → knee → ankle) that fold up
inside the wheel wells. The scene drives on its wheels, stops at a rock field,
swings the wheel modules outward, unfolds the legs, stands up, and walks a gait
in place while the cargo deck counter-rotates to stay roughly level.

Parts are each built in their OWN local frame; the runtime assembler
(yc_assemble) places them from the connector/mate plan, and every joint
origin/axis is derived with yc_pt / yc_dir so slider edits keep the rig joined.
Units mm, Z up.
"""
import math
import cadquery as cq

# ── YoCAD assembler (injected — do not edit) ──
def _yc_norm(v):
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n < 1e-12:
        return [0.0, 0.0, 1.0]
    return [v[0] / n, v[1] / n, v[2] / n]


def _yc_cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _yc_frame(c):
    """Connector dict -> [R(3x3 rows, columns = x,y,z basis), t]."""
    z = _yc_norm(c["z"])
    x = _yc_norm(c["x"])
    y = _yc_norm(_yc_cross(z, x))
    x = _yc_cross(y, z)
    R = [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]]
    o = c["origin"]
    return [R, [float(o[0]), float(o[1]), float(o[2])]]


def _yc_mul(A, B):
    Ra = A[0]
    Rb = B[0]
    R = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    i = 0
    while i < 3:
        j = 0
        while j < 3:
            R[i][j] = Ra[i][0] * Rb[0][j] + Ra[i][1] * Rb[1][j] + Ra[i][2] * Rb[2][j]
            j = j + 1
        i = i + 1
    tb = B[1]
    ta = A[1]
    t = [Ra[0][0] * tb[0] + Ra[0][1] * tb[1] + Ra[0][2] * tb[2] + ta[0],
         Ra[1][0] * tb[0] + Ra[1][1] * tb[1] + Ra[1][2] * tb[2] + ta[1],
         Ra[2][0] * tb[0] + Ra[2][1] * tb[1] + Ra[2][2] * tb[2] + ta[2]]
    return [R, t]


def _yc_inv(A):
    R = A[0]
    t = A[1]
    Rt = [[R[0][0], R[1][0], R[2][0]], [R[0][1], R[1][1], R[2][1]], [R[0][2], R[1][2], R[2][2]]]
    ti = [-(Rt[0][0] * t[0] + Rt[0][1] * t[1] + Rt[0][2] * t[2]),
          -(Rt[1][0] * t[0] + Rt[1][1] * t[1] + Rt[1][2] * t[2]),
          -(Rt[2][0] * t[0] + Rt[2][1] * t[1] + Rt[2][2] * t[2])]
    return [Rt, ti]


def _yc_axis_angle(R):
    tr = R[0][0] + R[1][1] + R[2][2]
    c = (tr - 1.0) / 2.0
    if c > 1.0:
        c = 1.0
    if c < -1.0:
        c = -1.0
    ang = math.acos(c)
    if ang < 1e-9:
        return [[0.0, 0.0, 1.0], 0.0]
    if ang > math.pi - 1e-6:
        xx = math.sqrt(max((R[0][0] + 1.0) / 2.0, 0.0))
        yy = math.sqrt(max((R[1][1] + 1.0) / 2.0, 0.0))
        zz = math.sqrt(max((R[2][2] + 1.0) / 2.0, 0.0))
        if xx >= yy and xx >= zz and xx > 1e-9:
            if R[0][1] < 0.0:
                yy = -yy
            if R[0][2] < 0.0:
                zz = -zz
            return [_yc_norm([xx, yy, zz]), 180.0]
        if yy >= zz and yy > 1e-9:
            if R[0][1] < 0.0:
                xx = -xx
            if R[1][2] < 0.0:
                zz = -zz
            return [_yc_norm([xx, yy, zz]), 180.0]
        if R[0][2] < 0.0:
            xx = -xx
        if R[1][2] < 0.0:
            yy = -yy
        return [_yc_norm([xx, yy, zz]), 180.0]
    s = 2.0 * math.sin(ang)
    ax = [(R[2][1] - R[1][2]) / s, (R[0][2] - R[2][0]) / s, (R[1][0] - R[0][1]) / s]
    return [_yc_norm(ax), ang * 180.0 / math.pi]


def yc_assemble(parts, connectors, mates, root):
    """Place parts by aligning mated connector frames. Returns [placed, xf]:
    placed = dict part -> transformed solid; xf = dict part -> world [R, t]."""
    flip = [[[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]], [0.0, 0.0, 0.0]]
    ident = [[[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.0, 0.0, 0.0]]
    xf = {}
    xf[root] = ident
    queue = [root]
    while len(queue) > 0:
        cur = queue.pop(0)
        k = 0
        while k < len(mates):
            a_id = mates[k][0]
            b_id = mates[k][1]
            ca = connectors[a_id]
            cb = connectors[b_id]
            pa = ca["part"]
            pb = cb["part"]
            if pa == cur and (pb in xf) == False:
                w = _yc_mul(_yc_mul(_yc_mul(xf[pa], _yc_frame(ca)), flip), _yc_inv(_yc_frame(cb)))
                xf[pb] = w
                queue.append(pb)
            if pb == cur and (pa in xf) == False:
                w = _yc_mul(_yc_mul(_yc_mul(xf[pb], _yc_frame(cb)), flip), _yc_inv(_yc_frame(ca)))
                xf[pa] = w
                queue.append(pa)
            k = k + 1
    placed = {}
    for name in parts:
        solid = parts[name]
        if name in xf:
            aa = _yc_axis_angle(xf[name][0])
            if aa[1] != 0.0:
                solid = solid.rotate((0, 0, 0), (aa[0][0], aa[0][1], aa[0][2]), aa[1])
            t = xf[name][1]
            solid = solid.translate((t[0], t[1], t[2]))
        placed[name] = solid
    return [placed, xf]


def yc_pt(xf, part, p):
    """Part-local point -> world (for JOINTS origins)."""
    A = xf[part]
    R = A[0]
    t = A[1]
    return [R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2] + t[0],
            R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2] + t[1],
            R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2] + t[2]]


def yc_dir(xf, part, v):
    """Part-local direction -> world (for JOINTS axes)."""
    R = xf[part][0]
    return _yc_norm([R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
                     R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
                     R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]])
# ── end YoCAD assembler ──


# ── Parameters (live sliders) ────────────────────────────────────────────────
CHASSIS_W = 1400.0     # central body width (X) [1000:1800]
CHASSIS_D = 1000.0     # central body depth (Y) [700:1300]
CHASSIS_H = 600.0      # central body height (Z) [420:800]
WHEEL_DIA = 520.0      # large wheel diameter [380:680]
WHEEL_W = 260.0        # wheel width [160:340]
THIGH_LEN = 350.0      # leg upper segment length [260:480]
SHIN_LEN = 350.0       # leg lower segment length [260:480]
DECK_W = 900.0         # cargo deck width [600:1150]
HALF_TRACK = 600.0     # wheel X offset from center [460:740]
HALF_BASE = 700.0      # wheel Y offset from center [520:860]
DEPLOY_MAX = 70.0      # wheel-module outward swing [30:90]

# ── Fixed proportions (derived / non-slider) ─────────────────────────────────
DECK_D = 700.0                     # cargo deck depth
DECK_H = 120.0                     # cargo deck thickness
MODULE_HIP_DZ = -100.0             # hip pivot height in wheel-module frame
MODULE_SPIN_DZ = -250.0            # wheel axle height in wheel-module frame
DEPLOY_Z = -CHASSIS_H / 4.0        # deploy pivot height on the chassis
DECK_MOUNT_Z = CHASSIS_H / 2.0     # cargo mount at rear top edge
# ride height so the wheels sit on the ground plane (z = 0) at rest
CHZ = -(DEPLOY_Z + MODULE_SPIN_DZ) + WHEEL_DIA / 2.0

BODY_W = 2.0 * HALF_TRACK
BODY_H = CHASSIS_H * 0.8
POD_R = 118.0
POD_LEN = 210.0

CORNERS = [("fl", -1, 1), ("fr", 1, 1), ("rl", -1, -1), ("rr", 1, -1)]

# ── Colors ───────────────────────────────────────────────────────────────────
COL_CHASSIS = cq.Color(0.90, 0.90, 0.92)
COL_DECK = cq.Color(0.85, 0.45, 0.10)
COL_MODULE = cq.Color(0.35, 0.38, 0.42)
COL_WHEEL = cq.Color(0.12, 0.12, 0.14)
COL_THIGH = cq.Color(0.80, 0.80, 0.82)
COL_SHIN = cq.Color(0.55, 0.57, 0.60)
COL_FOOT = cq.Color(0.20, 0.20, 0.22)


def _fillet(solid, sel, r):
    try:
        return solid.edges(sel).fillet(r)
    except Exception:
        return solid


# ── Chassis: pressurized body + fenders + deploy pods + rear deck mount ──────
def make_chassis():
    body = cq.Workplane("XY").box(BODY_W, CHASSIS_D, BODY_H)
    body = _fillet(body, "|Z", 110.0)

    # pressurized dome / cockpit hump on top
    dome = (cq.Workplane("XY").box(BODY_W * 0.55, CHASSIS_D * 0.6, BODY_H * 0.5)
            .translate((0.0, CHASSIS_D * 0.12, BODY_H * 0.42)))
    dome = _fillet(dome, "|Z", 120.0)
    dome = _fillet(dome, "#Z", 70.0)
    body = body.union(dome)

    # forward sensor bar
    body = body.union(
        cq.Workplane("XY").box(BODY_W * 0.5, 60.0, 60.0)
        .translate((0.0, CHASSIS_D / 2.0 - 30.0, BODY_H * 0.32)))
    # comms mast
    body = body.union(
        cq.Workplane("XY").circle(14.0).extrude(150.0)
        .translate((BODY_W * 0.28, CHASSIS_D * 0.18, BODY_H * 0.5)))

    for tag, sx, sy in CORNERS:
        cx = sx * HALF_TRACK
        cy = sy * HALF_BASE
        # fender arm reaching out to the wheel-module pivot
        arm_len = abs(cy) - CHASSIS_D / 2.0 + 180.0
        arm = (cq.Workplane("XY").box(220.0, arm_len, 260.0)
               .translate((cx, sy * (CHASSIS_D / 2.0 + arm_len / 2.0 - 120.0),
                           DEPLOY_Z)))
        arm = _fillet(arm, "|Z", 40.0)
        body = body.union(arm)
        # deploy-pivot barrel (axis along Y) — a readable hinge boss
        pod = (cq.Workplane("XZ").circle(POD_R).extrude(POD_LEN)
               .translate((cx, cy - POD_LEN / 2.0, DEPLOY_Z)))
        body = body.union(pod)

    # rear deck mount block so the cargo pivot sits ON the chassis
    body = body.union(
        cq.Workplane("XY").box(DECK_W * 0.7, 120.0, 140.0)
        .translate((0.0, -CHASSIS_D / 2.0 + 30.0, DECK_MOUNT_Z - 70.0)))

    return body.translate((0.0, 0.0, CHZ))


# ── Cargo deck: built HORIZONTAL in its own frame (X = width, -Y = depth,
#    Z = thickness/thin); the plan's mate keeps it upright, hinged along the
#    chassis rear-top edge so it reads as a flat payload bed, not a plate. ─────
def make_cargo_deck():
    # flat slab: thin in Z, extends rearward in -Y from the front hinge edge
    deck = (cq.Workplane("XY").box(DECK_W, DECK_D, DECK_H)
            .translate((0.0, -DECK_D / 2.0, -DECK_H / 2.0)))
    # tail wall rising at the far rear edge
    deck = deck.union(
        cq.Workplane("XY").box(DECK_W, 70.0, DECK_H + 150.0)
        .translate((0.0, -DECK_D + 35.0, 75.0 - DECK_H / 2.0)))
    # two side rails running along the deck length
    for s in (-1.0, 1.0):
        deck = deck.union(
            cq.Workplane("XY").box(60.0, DECK_D, DECK_H + 90.0)
            .translate((s * (DECK_W / 2.0 - 30.0), -DECK_D / 2.0,
                        45.0 - DECK_H / 2.0)))
    return deck


# ── Wheel module: hangs down from the deploy pivot; symmetric so no handedness.
#    Local Z up; wm_deploy at 0, wm_hip at -100, wm_spin (axle) at -250. ───────
def make_wheel_module():
    housing = (cq.Workplane("XY").box(300.0, 360.0, 380.0)
               .translate((0.0, 0.0, -140.0)))
    housing = _fillet(housing, "|Y", 46.0)
    # deploy knuckle barrel (axis local Y) at z = 0
    housing = housing.union(
        cq.Workplane("XZ").circle(112.0).extrude(300.0)
        .translate((0.0, -150.0, 0.0)))
    # hollow wheel well on the -Y face where the folded leg tucks
    housing = housing.cut(
        cq.Workplane("XY").box(220.0, 200.0, 300.0)
        .translate((0.0, -150.0, -160.0)))
    # axle hub boss (axis local X) at the wheel line
    housing = housing.union(
        cq.Workplane("YZ").circle(74.0).extrude(WHEEL_W + 60.0)
        .translate((-(WHEEL_W + 60.0) / 2.0, 0.0, MODULE_SPIN_DZ)))
    # hip pivot boss (axis local X)
    housing = housing.union(
        cq.Workplane("YZ").circle(58.0).extrude(220.0)
        .translate((-110.0, 0.0, MODULE_HIP_DZ)))
    return housing


# ── Wheel: spins about local X, symmetric ────────────────────────────────────
def make_wheel():
    r = WHEEL_DIA / 2.0
    w = WHEEL_W
    tire = (cq.Workplane("YZ").circle(r).extrude(w).translate((-w / 2.0, 0.0, 0.0)))
    tire = tire.cut(cq.Workplane("YZ").circle(r * 0.82).extrude(w * 0.5)
                    .translate((-w * 0.25, 0.0, 0.0)))
    tire = tire.cut(cq.Workplane("YZ").circle(r * 0.82).extrude(w * 0.5)
                    .translate((-w * 0.25, 0.0, 0.0)).translate((w * 0.5, 0, 0)))
    hub = (cq.Workplane("YZ").circle(r * 0.34).extrude(w * 0.7)
           .translate((-w * 0.35, 0.0, 0.0)))
    wheel = tire.union(hub)
    # spokes
    for i in range(5):
        a = math.pi * i / 2.5
        wheel = wheel.union(
            cq.Workplane("XY").box(w * 0.5, 44.0, r * 1.5)
            .rotate((0, 0, 0), (1, 0, 0), math.degrees(a)))
    # tread lugs
    n = 10
    for i in range(n):
        a = 2.0 * math.pi * i / n
        lug = (cq.Workplane("XY").box(w - 34.0, 48.0, 46.0)
               .translate((0.0, 0.0, r - 6.0))
               .rotate((0, 0, 0), (1, 0, 0), math.degrees(a)))
        wheel = wheel.union(lug)
    return wheel


def _limb(cross, length):
    seg = (cq.Workplane("XY").box(cross, cross, length)
           .translate((0.0, 0.0, -length / 2.0)))
    # joint knuckles (axis local X) at both ends make the pivots readable
    for z in (0.0, -length):
        seg = seg.union(
            cq.Workplane("YZ").circle(cross * 0.55).extrude(cross + 40.0)
            .translate((-(cross + 40.0) / 2.0, 0.0, z)))
    return seg


def make_thigh():
    return _limb(120.0, THIGH_LEN)


def make_shin():
    return _limb(100.0, SHIN_LEN)


def make_foot():
    ankle = (cq.Workplane("YZ").circle(58.0).extrude(150.0)
             .translate((-75.0, 0.0, 0.0)))
    pad = (cq.Workplane("XY").box(220.0, 220.0, 90.0)
           .translate((0.0, 0.0, -45.0)))
    pad = _fillet(pad, "|Z", 30.0)
    # sole grip ridges
    for gy in (-70.0, 0.0, 70.0):
        pad = pad.union(cq.Workplane("XY").box(200.0, 30.0, 24.0)
                        .translate((0.0, gy, -90.0)))
    return ankle.union(pad)


# ── Build every part in its local frame ──────────────────────────────────────
PARTS = {"chassis": make_chassis(), "cargo_deck": make_cargo_deck()}
for tag, sx, sy in CORNERS:
    PARTS["wheel_module_" + tag] = make_wheel_module()
    PARTS["wheel_" + tag] = make_wheel()
    PARTS["thigh_" + tag] = make_thigh()
    PARTS["shin_" + tag] = make_shin()
    PARTS["foot_" + tag] = make_foot()

# ── Connector frames (from the plan; z_axis→"z", x_axis→"x") ─────────────────
CONNECTORS = {
    "chassis_cargo": {"part": "chassis",
                      "origin": [0.0, -CHASSIS_D / 2.0, DECK_MOUNT_Z + CHZ],
                      "z": [0, 0, 1], "x": [1, 0, 0]},
    "deck_root": {"part": "cargo_deck", "origin": [0.0, 0.0, 0.0],
                  "z": [0, 0, -1], "x": [1, 0, 0]},
}
for tag, sx, sy in CORNERS:
    wm = "wheel_module_" + tag
    CONNECTORS["chassis_deploy_" + tag] = {
        "part": "chassis",
        "origin": [sx * HALF_TRACK, sy * HALF_BASE, DEPLOY_Z + CHZ],
        "z": [0, 1, 0], "x": [0, 0, 1]}
    CONNECTORS["wm_deploy_" + tag] = {"part": wm, "origin": [0.0, 0.0, 0.0],
                                      "z": [0, 1, 0], "x": [0, 0, 1]}
    CONNECTORS["wm_spin_" + tag] = {"part": wm, "origin": [0.0, 0.0, MODULE_SPIN_DZ],
                                    "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["wm_hip_" + tag] = {"part": wm, "origin": [0.0, 0.0, MODULE_HIP_DZ],
                                   "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["wheel_spin_" + tag] = {"part": "wheel_" + tag, "origin": [0.0, 0.0, 0.0],
                                       "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["thigh_hip_" + tag] = {"part": "thigh_" + tag, "origin": [0.0, 0.0, 0.0],
                                      "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["thigh_knee_" + tag] = {"part": "thigh_" + tag,
                                       "origin": [0.0, 0.0, -THIGH_LEN],
                                       "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["shin_knee_" + tag] = {"part": "shin_" + tag, "origin": [0.0, 0.0, 0.0],
                                      "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["shin_ankle_" + tag] = {"part": "shin_" + tag,
                                       "origin": [0.0, 0.0, -SHIN_LEN],
                                       "z": [1, 0, 0], "x": [0, 0, 1]}
    CONNECTORS["foot_ankle_" + tag] = {"part": "foot_" + tag, "origin": [0.0, 0.0, 0.0],
                                       "z": [1, 0, 0], "x": [0, 0, 1]}

# ── Mates (from the plan) ────────────────────────────────────────────────────
MATES = [["chassis_cargo", "deck_root"]]
for tag, sx, sy in CORNERS:
    MATES.append(["chassis_deploy_" + tag, "wm_deploy_" + tag])
    MATES.append(["wm_spin_" + tag, "wheel_spin_" + tag])
    MATES.append(["wm_hip_" + tag, "thigh_hip_" + tag])
    MATES.append(["thigh_knee_" + tag, "shin_knee_" + tag])
    MATES.append(["shin_ankle_" + tag, "foot_ankle_" + tag])

# ── Deterministic placement ──────────────────────────────────────────────────
_asm = yc_assemble(PARTS, CONNECTORS, MATES, "chassis")
PLACED = _asm[0]
XF = _asm[1]

result = (cq.Assembly()
          .add(PLACED["chassis"], name="chassis", color=COL_CHASSIS)
          .add(PLACED["cargo_deck"], name="cargo_deck", color=COL_DECK))
for tag, sx, sy in CORNERS:
    result = result.add(PLACED["wheel_module_" + tag], name="wheel_module_" + tag, color=COL_MODULE)
    result = result.add(PLACED["wheel_" + tag], name="wheel_" + tag, color=COL_WHEEL)
    result = result.add(PLACED["thigh_" + tag], name="thigh_" + tag, color=COL_THIGH)
    result = result.add(PLACED["shin_" + tag], name="shin_" + tag, color=COL_SHIN)
    result = result.add(PLACED["foot_" + tag], name="foot_" + tag, color=COL_FOOT)

# ── Joints (origins/axes derived from the assembled transforms) ──────────────
JOINTS = [
    {"name": "cargo_level", "type": "revolute", "parts": ["cargo_deck", "chassis"],
     "origin": yc_pt(XF, "chassis", [0.0, -CHASSIS_D / 2.0, DECK_MOUNT_Z + CHZ]),
     "axis": yc_dir(XF, "chassis", [1, 0, 0]), "limits": [-20.0, 20.0]},
]
for tag, sx, sy in CORNERS:
    wm = "wheel_module_" + tag
    JOINTS.append({"name": "deploy_" + tag, "type": "revolute",
                   "parts": [wm, "chassis"],
                   "origin": yc_pt(XF, "chassis",
                                   [sx * HALF_TRACK, sy * HALF_BASE, DEPLOY_Z + CHZ]),
                   "axis": yc_dir(XF, "chassis", [0, 1, 0]), "limits": [0.0, 90.0]})
    JOINTS.append({"name": "spin_" + tag, "type": "revolute",
                   "parts": ["wheel_" + tag, wm],
                   "origin": yc_pt(XF, wm, [0.0, 0.0, MODULE_SPIN_DZ]),
                   "axis": yc_dir(XF, wm, [1, 0, 0]), "limits": [0.0, 360.0]})
    JOINTS.append({"name": "hip_" + tag, "type": "revolute",
                   "parts": ["thigh_" + tag, wm],
                   "origin": yc_pt(XF, wm, [0.0, 0.0, MODULE_HIP_DZ]),
                   "axis": yc_dir(XF, wm, [1, 0, 0]), "limits": [-45.0, 45.0]})
    JOINTS.append({"name": "knee_" + tag, "type": "revolute",
                   "parts": ["shin_" + tag, "thigh_" + tag],
                   "origin": yc_pt(XF, "thigh_" + tag, [0.0, 0.0, -THIGH_LEN]),
                   "axis": yc_dir(XF, "thigh_" + tag, [1, 0, 0]), "limits": [0.0, 120.0]})
    JOINTS.append({"name": "ankle_" + tag, "type": "revolute",
                   "parts": ["foot_" + tag, "shin_" + tag],
                   "origin": yc_pt(XF, "shin_" + tag, [0.0, 0.0, -SHIN_LEN]),
                   "axis": yc_dir(XF, "shin_" + tag, [1, 0, 0]), "limits": [-45.0, 45.0]})

LINKS = []

# ── Animation ────────────────────────────────────────────────────────────────
# Sequence: drive (0-3) → stop (3-4) → swing modules out (4-6) → unfold & stand
# (6-8) → walk gait in place with the deck actively leveling (8-12) → fold back
# to the driving pose (12-15). Diagonal legs (fl,rr) and (fr,rl) step out of phase.
HIP_TUCK, KNEE_TUCK, ANK_TUCK = 28.0, 112.0, -28.0   # folded up in the wheel well
HIP_DOWN, KNEE_DOWN, ANK_DOWN = -6.0, 74.0, -34.0    # legs planted as stabilizers

TRACKS = []

# Wheels spin about their own transverse axle while driving in, and again after
# the stabilizers retract as the rover drives off. The wheel modules NEVER swing
# out, so all four wheels stay on the ground the entire time.
for tag, sx, sy in CORNERS:
    TRACKS.append({"joint": "spin_" + tag, "unit": "deg",
                   "keys": [[0.0, 0.0], [3.0, 360.0], [11.0, 360.0],
                            [15.0, 720.0]]})

# All four legs deploy straight DOWN as outrigger stabilizers with IDENTICAL
# angles, so the stance is fully symmetric and the feet drop toward the ground
# beside the wheels (no lateral splay, no lifting the wheels). Then they fold
# back into the wheel wells. The hip axis is transverse, so the legs move only
# fore/aft-and-down within the ground plane — never up into the air.
for tag, sx, sy in CORNERS:
    TRACKS.append({"joint": "hip_" + tag, "unit": "deg",
                   "keys": [[0.0, HIP_TUCK], [4.0, HIP_TUCK], [6.0, HIP_DOWN],
                            [9.0, HIP_DOWN], [11.0, HIP_TUCK], [15.0, HIP_TUCK]]})
    TRACKS.append({"joint": "knee_" + tag, "unit": "deg",
                   "keys": [[0.0, KNEE_TUCK], [4.0, KNEE_TUCK], [6.0, KNEE_DOWN],
                            [9.0, KNEE_DOWN], [11.0, KNEE_TUCK], [15.0, KNEE_TUCK]]})
    TRACKS.append({"joint": "ankle_" + tag, "unit": "deg",
                   "keys": [[0.0, ANK_TUCK], [4.0, ANK_TUCK], [6.0, ANK_DOWN],
                            [9.0, ANK_DOWN], [11.0, ANK_TUCK], [15.0, ANK_TUCK]]})

# cargo deck trims gently to stay level while the rover is parked on its legs
TRACKS.append({"joint": "cargo_level", "unit": "deg",
               "keys": [[0.0, 0.0], [6.0, 0.0], [7.5, 5.0], [9.0, -5.0],
                        [10.5, 0.0], [15.0, 0.0]]})

ANIMATION = {"fps": 24, "duration": 15.0, "tracks": TRACKS}