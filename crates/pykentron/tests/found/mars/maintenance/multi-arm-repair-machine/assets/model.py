"""Compact Mars maintenance robot — a wheeled mobile base (six rocker wheels
folded into the chassis) carrying four deployable outrigger stabilizers, a
rotating upper turret, and three articulated repair arms. Each arm is a serial
shoulder→elbow→wrist chain ending in an interchangeable tool: a parallel
gripper, an inspection camera, and a powered wrench with a spinning socket.

The scene is built entirely in local part frames and placed by the deterministic
connector assembler (yc_assemble), exactly from the supplied kinematic plan;
every JOINTS origin/axis is derived with yc_pt/yc_dir so slider edits keep the
whole robot attached. The timeline plays the concept mission: park beside a
damaged equipment panel, deploy the stabilizers, reach in with the gripper arm
to open and hold the panel, swing the camera arm in to inspect the interior, and
finally drive the wrench arm down to spin-tighten one fastener before stowing.

Concept only — not a validated repair system.
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


# ── Parameters ───────────────────────────────────────────────────────────────
BASE_W = 600.0          # chassis width (X) [400:800]
BASE_D = 600.0          # chassis depth (Y) [400:800]
BASE_H = 260.0          # chassis body height [180:340]
WHEEL_R = 100.0         # rover wheel radius [70:140]
TURRET_R = 130.0        # turret radius [90:170]
TURRET_H = 200.0        # turret height [120:280]
UPPER_ARM_LEN = 300.0   # arm shoulder→elbow length [200:400]
FOREARM_LEN = 250.0     # arm elbow→wrist length [160:340]
STAB_LEN = 320.0        # stabilizer leg length [220:440]
STAB_DEPLOY = 75.0      # stabilizer deploy angle, deg [40:90]

# ── fixed proportions / derived layout ──────────────────────────────────────
WRIST_LEN = 90.0                      # wrist link length (tool standoff)
WHEEL_Z = 60.0                        # wheel axle height on the chassis side
GROUND = WHEEL_Z - WHEEL_R            # ground plane (wheel contact)
SH_Z = TURRET_H * 0.55               # shoulder height on the turret
SH_R = TURRET_R - 10.0                # shoulder ring radius on the turret
STAB_XO = BASE_W / 2.0 - 50.0         # stabilizer mount X offset
STAB_YO = BASE_D / 2.0 - 50.0         # stabilizer mount Y offset
STAB_Z = 90.0                         # stabilizer mount height on the base side
COS120 = -0.5
SIN120 = 0.8660254

# shoulder ring positions (turret local frame, 120° apart)
SH1 = (SH_R, 0.0, SH_Z)
SH2 = (SH_R * COS120, SH_R * SIN120, SH_Z)
SH3 = (SH_R * COS120, -SH_R * SIN120, SH_Z)
# radial "out" direction of each shoulder (turret local)
RAD1 = (1.0, 0.0, 0.0)
RAD2 = (COS120, SIN120, 0.0)
RAD3 = (COS120, -SIN120, 0.0)

# ── colors (from the plan) ───────────────────────────────────────────────────
COLORS = {
    "mobile_base": cq.Color(0.55, 0.55, 0.58),
    "turret": cq.Color(0.80, 0.80, 0.83),
    "stabilizer_fl": cq.Color(0.85, 0.55, 0.15),
    "stabilizer_fr": cq.Color(0.85, 0.55, 0.15),
    "stabilizer_rl": cq.Color(0.85, 0.55, 0.15),
    "stabilizer_rr": cq.Color(0.85, 0.55, 0.15),
    "upper_arm_1": cq.Color(0.75, 0.76, 0.78),
    "forearm_1": cq.Color(0.70, 0.71, 0.73),
    "wrist_1": cq.Color(0.60, 0.61, 0.63),
    "gripper_base_1": cq.Color(0.20, 0.55, 0.85),
    "gripper_finger_1": cq.Color(0.20, 0.55, 0.85),
    "upper_arm_2": cq.Color(0.75, 0.76, 0.78),
    "forearm_2": cq.Color(0.70, 0.71, 0.73),
    "wrist_2": cq.Color(0.60, 0.61, 0.63),
    "inspection_camera": cq.Color(0.15, 0.65, 0.45),
    "upper_arm_3": cq.Color(0.75, 0.76, 0.78),
    "forearm_3": cq.Color(0.70, 0.71, 0.73),
    "wrist_3": cq.Color(0.60, 0.61, 0.63),
    "powered_wrench": cq.Color(0.85, 0.25, 0.25),
}


# ── local-frame part builders ────────────────────────────────────────────────
def _box(w, d, h, z0=0.0):
    """Box spanning z0..z0+h, centered in X and Y."""
    return cq.Workplane("XY").box(w, d, h).translate((0.0, 0.0, z0 + h / 2.0))


def _safe_fillet(solid, r):
    # Fillets omitted on purpose: across a ~20-part articulated assembly the
    # extra curved edges multiply the tessellation/collision face count enough
    # to time out the native geometry check. Chamfer-free prisms read fine here.
    return solid


def make_base():
    """Wheeled chassis, root part. Body spans z 0..BASE_H; six wheels folded
    into the ±Y sides sit the rover on the ground at z=GROUND."""
    body = _box(BASE_W, BASE_D, BASE_H)
    body = _safe_fillet(body, 24.0)
    # sensor deck lip around the top edge
    body = body.union(_box(BASE_W - 40.0, BASE_D - 40.0, 18.0, BASE_H - 18.0))
    # six rocker wheels (3 per side), axis along Y
    for xw in (-BASE_W / 2.0 + 120.0, 0.0, BASE_W / 2.0 - 120.0):
        for yw in (BASE_D / 2.0, -BASE_D / 2.0):
            wheel = (cq.Workplane("XY").circle(WHEEL_R).extrude(70.0)
                     .translate((0.0, 0.0, -35.0))
                     .rotate((0, 0, 0), (1, 0, 0), 90.0)
                     .translate((xw, yw, WHEEL_Z)))
            body = body.union(wheel)
    return body


def make_turret():
    """Rotating turret ring; local origin at its base (mates onto the chassis
    top). Shoulder bosses ring the upper edge 120° apart."""
    t = cq.Workplane("XY").circle(TURRET_R).extrude(TURRET_H)
    t = t.union(cq.Workplane("XY").circle(TURRET_R + 6.0).extrude(16.0))
    t = t.union(cq.Workplane("XY").circle(TURRET_R - 18.0).extrude(24.0)
                .translate((0.0, 0.0, TURRET_H - 4.0)))  # sink 4 mm -> real overlap
    for pos, rad in ((SH1, RAD1), (SH2, RAD2), (SH3, RAD3)):
        # box boss (prism-into-cylinder union is far cheaper for the native core
        # than the old cylinder-into-cylinder fuse, and reads as a mount pad)
        boss = (cq.Workplane("XY").box(60.0, 60.0, 60.0)
                .rotate((0, 0, 0), (0, 0, 1), math.degrees(math.atan2(rad[1], rad[0])))
                .translate(pos))
        t = t.union(boss)
    return t


def make_stabilizer():
    """Deployable outrigger leg; local origin at the pivot, leg extends +Z
    (stowed upright) with a wide foot pad at the far end."""
    leg = _box(60.0, 60.0, STAB_LEN)
    leg = _safe_fillet(leg, 10.0)
    leg = leg.union(_box(110.0, 110.0, 22.0, STAB_LEN - 22.0))  # foot pad
    return leg


def make_upper_arm():
    a = _box(70.0, 70.0, UPPER_ARM_LEN)
    # boxy elbow boss: a plain prism overlap keeps the joint area readable
    # without the curved cylinder-box fuse that inflates the geometry check
    a = a.union(_box(80.0, 76.0, 60.0, UPPER_ARM_LEN - 30.0))
    return a


def make_forearm():
    a = _box(60.0, 60.0, FOREARM_LEN)
    a = a.union(_box(70.0, 66.0, 52.0, FOREARM_LEN - 26.0))  # boxy wrist boss
    return a


def make_wrist():
    w = _box(55.0, 55.0, WRIST_LEN)
    w = _safe_fillet(w, 8.0)
    return w


def make_gripper_base():
    """Parallel-gripper body; a fixed jaw on -Y and a mount for the moving
    finger on +Y (finger connector at local [0,25,70])."""
    g = _box(70.0, 70.0, 90.0)
    g = _safe_fillet(g, 8.0)
    g = g.union(_box(20.0, 16.0, 64.0).translate((0.0, -27.0, 88.0)))  # fixed jaw (2 mm overlap)
    return g


def make_gripper_finger():
    return _safe_fillet(_box(20.0, 50.0, 80.0), 6.0)


def make_camera():
    c = _box(70.0, 70.0, 80.0)
    c = c.union(cq.Workplane("XY").circle(22.0).extrude(28.0).translate((0.0, 0.0, 76.0)))
    return c


def make_wrench():
    w = _box(70.0, 70.0, 80.0)
    w = _safe_fillet(w, 8.0)
    w = w.union(cq.Workplane("XY").circle(24.0).extrude(24.0).translate((0.0, 0.0, 76.0)))
    w = w.union(cq.Workplane("XY").polygon(6, 42.0).extrude(14.0).translate((0.0, 0.0, 96.0)))  # socket
    return w


# ── PARTS (local frames) ─────────────────────────────────────────────────────
PARTS = {
    "mobile_base": make_base(),
    "turret": make_turret(),
    "stabilizer_fl": make_stabilizer(),
    "stabilizer_fr": make_stabilizer(),
    "stabilizer_rl": make_stabilizer(),
    "stabilizer_rr": make_stabilizer(),
    "upper_arm_1": make_upper_arm(),
    "forearm_1": make_forearm(),
    "wrist_1": make_wrist(),
    "gripper_base_1": make_gripper_base(),
    "gripper_finger_1": make_gripper_finger(),
    "upper_arm_2": make_upper_arm(),
    "forearm_2": make_forearm(),
    "wrist_2": make_wrist(),
    "inspection_camera": make_camera(),
    "upper_arm_3": make_upper_arm(),
    "forearm_3": make_forearm(),
    "wrist_3": make_wrist(),
    "powered_wrench": make_wrench(),
}

# ── CONNECTORS (from the plan; origins expressed via named parameters) ───────
CONNECTORS = {
    "c_base_turret":  {"part": "mobile_base", "origin": [0.0, 0.0, BASE_H], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_base_stab_fl": {"part": "mobile_base", "origin": [STAB_XO, STAB_YO, STAB_Z], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_base_stab_fr": {"part": "mobile_base", "origin": [STAB_XO, -STAB_YO, STAB_Z], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_base_stab_rl": {"part": "mobile_base", "origin": [-STAB_XO, STAB_YO, STAB_Z], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_base_stab_rr": {"part": "mobile_base", "origin": [-STAB_XO, -STAB_YO, STAB_Z], "z": [0, 0, 1], "x": [1, 0, 0]},

    "c_turret_base": {"part": "turret", "origin": [0.0, 0.0, 0.0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_turret_sh1":  {"part": "turret", "origin": [SH1[0], SH1[1], SH1[2]], "z": [RAD1[0], RAD1[1], 0.0], "x": [0, 0, 1]},
    "c_turret_sh2":  {"part": "turret", "origin": [SH2[0], SH2[1], SH2[2]], "z": [RAD2[0], RAD2[1], 0.0], "x": [0, 0, 1]},
    "c_turret_sh3":  {"part": "turret", "origin": [SH3[0], SH3[1], SH3[2]], "z": [RAD3[0], RAD3[1], 0.0], "x": [0, 0, 1]},

    "c_stab_fl_root": {"part": "stabilizer_fl", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_stab_fr_root": {"part": "stabilizer_fr", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_stab_rl_root": {"part": "stabilizer_rl", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_stab_rr_root": {"part": "stabilizer_rr", "origin": [0, 0, 0], "z": [0, 0, 1], "x": [1, 0, 0]},

    "c_ua1_root":  {"part": "upper_arm_1", "origin": [0, 0, 0], "z": [-1, 0, 0], "x": [0, 0, 1]},
    "c_ua1_elbow": {"part": "upper_arm_1", "origin": [0.0, 0.0, UPPER_ARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_fa1_root":  {"part": "forearm_1", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_fa1_wrist": {"part": "forearm_1", "origin": [0.0, 0.0, FOREARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wr1_root":  {"part": "wrist_1", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_wr1_tool":  {"part": "wrist_1", "origin": [0.0, 0.0, WRIST_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_gr1_root":   {"part": "gripper_base_1", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_gr1_finger": {"part": "gripper_base_1", "origin": [0, 25, 70], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_gf1_root":   {"part": "gripper_finger_1", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},

    "c_ua2_root":  {"part": "upper_arm_2", "origin": [0, 0, 0], "z": [-1, 0, 0], "x": [0, 0, 1]},
    "c_ua2_elbow": {"part": "upper_arm_2", "origin": [0.0, 0.0, UPPER_ARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_fa2_root":  {"part": "forearm_2", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_fa2_wrist": {"part": "forearm_2", "origin": [0.0, 0.0, FOREARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wr2_root":  {"part": "wrist_2", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_wr2_tool":  {"part": "wrist_2", "origin": [0.0, 0.0, WRIST_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_cam2_root": {"part": "inspection_camera", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},

    "c_ua3_root":  {"part": "upper_arm_3", "origin": [0, 0, 0], "z": [-1, 0, 0], "x": [0, 0, 1]},
    "c_ua3_elbow": {"part": "upper_arm_3", "origin": [0.0, 0.0, UPPER_ARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_fa3_root":  {"part": "forearm_3", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_fa3_wrist": {"part": "forearm_3", "origin": [0.0, 0.0, FOREARM_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wr3_root":  {"part": "wrist_3", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
    "c_wr3_tool":  {"part": "wrist_3", "origin": [0.0, 0.0, WRIST_LEN], "z": [0, 0, 1], "x": [1, 0, 0]},
    "c_wrench3_root": {"part": "powered_wrench", "origin": [0, 0, 0], "z": [0, 0, -1], "x": [1, 0, 0]},
}

# ── MATES (from the plan) ────────────────────────────────────────────────────
MATES = [
    ["c_base_turret", "c_turret_base"],
    ["c_base_stab_fl", "c_stab_fl_root"],
    ["c_base_stab_fr", "c_stab_fr_root"],
    ["c_base_stab_rl", "c_stab_rl_root"],
    ["c_base_stab_rr", "c_stab_rr_root"],
    ["c_turret_sh1", "c_ua1_root"],
    ["c_ua1_elbow", "c_fa1_root"],
    ["c_fa1_wrist", "c_wr1_root"],
    ["c_wr1_tool", "c_gr1_root"],
    ["c_gr1_finger", "c_gf1_root"],
    ["c_turret_sh2", "c_ua2_root"],
    ["c_ua2_elbow", "c_fa2_root"],
    ["c_fa2_wrist", "c_wr2_root"],
    ["c_wr2_tool", "c_cam2_root"],
    ["c_turret_sh3", "c_ua3_root"],
    ["c_ua3_elbow", "c_fa3_root"],
    ["c_fa3_wrist", "c_wr3_root"],
    ["c_wr3_tool", "c_wrench3_root"],
]

# ── deterministic placement ──────────────────────────────────────────────────
_asm = yc_assemble(PARTS, CONNECTORS, MATES, "mobile_base")
PLACED = _asm[0]
XF = _asm[1]

# ── result assembly ──────────────────────────────────────────────────────────
result = cq.Assembly()
for _n in PLACED:
    result = result.add(PLACED[_n], name=_n, color=COLORS[_n])

# static context prop: the damaged equipment panel the robot services (grounded,
# off the +X side where the gripper arm reaches). Not part of the kinematic chain.
PANEL_X = BASE_W / 2.0 + 240.0
PANEL_H = 440.0
panel = _box(70.0, 320.0, PANEL_H, GROUND).translate((PANEL_X, 0.0, 0.0))
# recessed access hatch on the face toward the robot
panel = panel.cut(_box(30.0, 220.0, 300.0, GROUND + 70.0)
                  .translate((PANEL_X - 30.0, 0.0, 0.0)))
result = result.add(panel, name="equipment_panel", color=cq.Color(0.42, 0.44, 0.48))

# ── joints (all origins/axes derived via yc_pt / yc_dir) ─────────────────────
JOINTS = [
    {"name": "turret_yaw", "type": "revolute", "parts": ["turret", "mobile_base"],
     "origin": yc_pt(XF, "turret", [0, 0, 0]), "axis": yc_dir(XF, "mobile_base", [0, 0, 1]),
     "limits": [-180.0, 180.0]},

    {"name": "stabilizer_fl_deploy", "type": "revolute", "parts": ["stabilizer_fl", "mobile_base"],
     "origin": yc_pt(XF, "stabilizer_fl", [0, 0, 0]), "axis": yc_dir(XF, "mobile_base", [-0.707, 0.707, 0]),
     "limits": [0.0, 90.0]},
    {"name": "stabilizer_fr_deploy", "type": "revolute", "parts": ["stabilizer_fr", "mobile_base"],
     "origin": yc_pt(XF, "stabilizer_fr", [0, 0, 0]), "axis": yc_dir(XF, "mobile_base", [0.707, 0.707, 0]),
     "limits": [0.0, 90.0]},
    {"name": "stabilizer_rl_deploy", "type": "revolute", "parts": ["stabilizer_rl", "mobile_base"],
     "origin": yc_pt(XF, "stabilizer_rl", [0, 0, 0]), "axis": yc_dir(XF, "mobile_base", [-0.707, -0.707, 0]),
     "limits": [0.0, 90.0]},
    {"name": "stabilizer_rr_deploy", "type": "revolute", "parts": ["stabilizer_rr", "mobile_base"],
     "origin": yc_pt(XF, "stabilizer_rr", [0, 0, 0]), "axis": yc_dir(XF, "mobile_base", [0.707, -0.707, 0]),
     "limits": [0.0, 90.0]},

    # arm 1 — gripper
    {"name": "shoulder_1", "type": "revolute", "parts": ["upper_arm_1", "turret"],
     "origin": yc_pt(XF, "upper_arm_1", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_1", [0, 1, 0]),
     "limits": [-95.0, 95.0]},
    {"name": "elbow_1", "type": "revolute", "parts": ["forearm_1", "upper_arm_1"],
     "origin": yc_pt(XF, "forearm_1", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_1", [0, 1, 0]),
     "limits": [0.0, 150.0]},
    {"name": "wrist_1", "type": "revolute", "parts": ["wrist_1", "forearm_1"],
     "origin": yc_pt(XF, "wrist_1", [0, 0, 0]), "axis": yc_dir(XF, "forearm_1", [1, 0, 0]),
     "limits": [-120.0, 120.0]},
    {"name": "tool_mount_1", "type": "fixed", "parts": ["gripper_base_1", "wrist_1"],
     "origin": yc_pt(XF, "gripper_base_1", [0, 0, 0]), "axis": yc_dir(XF, "wrist_1", [0, 0, 1]),
     "limits": [0.0, 0.0]},
    {"name": "gripper_1", "type": "revolute", "parts": ["gripper_finger_1", "gripper_base_1"],
     "origin": yc_pt(XF, "gripper_finger_1", [0, 0, 0]), "axis": yc_dir(XF, "gripper_base_1", [1, 0, 0]),
     "limits": [0.0, 40.0]},

    # arm 2 — inspection camera
    {"name": "shoulder_2", "type": "revolute", "parts": ["upper_arm_2", "turret"],
     "origin": yc_pt(XF, "upper_arm_2", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_2", [0, 1, 0]),
     "limits": [-95.0, 95.0]},
    {"name": "elbow_2", "type": "revolute", "parts": ["forearm_2", "upper_arm_2"],
     "origin": yc_pt(XF, "forearm_2", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_2", [0, 1, 0]),
     "limits": [0.0, 150.0]},
    {"name": "wrist_2", "type": "revolute", "parts": ["wrist_2", "forearm_2"],
     "origin": yc_pt(XF, "wrist_2", [0, 0, 0]), "axis": yc_dir(XF, "forearm_2", [1, 0, 0]),
     "limits": [-120.0, 120.0]},
    {"name": "tool_mount_2", "type": "fixed", "parts": ["inspection_camera", "wrist_2"],
     "origin": yc_pt(XF, "inspection_camera", [0, 0, 0]), "axis": yc_dir(XF, "wrist_2", [0, 0, 1]),
     "limits": [0.0, 0.0]},

    # arm 3 — powered wrench
    {"name": "shoulder_3", "type": "revolute", "parts": ["upper_arm_3", "turret"],
     "origin": yc_pt(XF, "upper_arm_3", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_3", [0, 1, 0]),
     "limits": [-95.0, 95.0]},
    {"name": "elbow_3", "type": "revolute", "parts": ["forearm_3", "upper_arm_3"],
     "origin": yc_pt(XF, "forearm_3", [0, 0, 0]), "axis": yc_dir(XF, "upper_arm_3", [0, 1, 0]),
     "limits": [0.0, 150.0]},
    {"name": "wrist_3", "type": "revolute", "parts": ["wrist_3", "forearm_3"],
     "origin": yc_pt(XF, "wrist_3", [0, 0, 0]), "axis": yc_dir(XF, "forearm_3", [1, 0, 0]),
     "limits": [-120.0, 120.0]},
    {"name": "wrench_spin_3", "type": "revolute", "parts": ["powered_wrench", "wrist_3"],
     "origin": yc_pt(XF, "powered_wrench", [0, 0, 0]), "axis": yc_dir(XF, "powered_wrench", [0, 0, 1]),
     "limits": [0.0, 360.0]},
]
LINKS = []

# ── animation: the concept repair mission (16 s) ─────────────────────────────
# Every excursion is anchored at rest before and after so nothing drifts across
# the gap to its next key. Values stay inside each joint's limits.
ANIMATION = {"fps": 24, "duration": 16.0, "tracks": [
    # 1) deploy the four stabilizers, hold through the mission, then stow
    {"joint": "stabilizer_fl_deploy", "unit": "deg",
     "keys": [[0, 0], [1.5, 0], [3.5, STAB_DEPLOY], [13.5, STAB_DEPLOY], [15.5, 0], [16, 0]]},
    {"joint": "stabilizer_fr_deploy", "unit": "deg",
     "keys": [[0, 0], [1.5, 0], [3.5, STAB_DEPLOY], [13.5, STAB_DEPLOY], [15.5, 0], [16, 0]]},
    {"joint": "stabilizer_rl_deploy", "unit": "deg",
     "keys": [[0, 0], [1.5, 0], [3.5, STAB_DEPLOY], [13.5, STAB_DEPLOY], [15.5, 0], [16, 0]]},
    {"joint": "stabilizer_rr_deploy", "unit": "deg",
     "keys": [[0, 0], [1.5, 0], [3.5, STAB_DEPLOY], [13.5, STAB_DEPLOY], [15.5, 0], [16, 0]]},

    # 2) turret rotates the gripper arm to face the panel, holds, returns
    {"joint": "turret_yaw", "unit": "deg",
     "keys": [[0, 0], [3.5, 0], [5, 40], [13.5, 40], [15.5, 0], [16, 0]]},

    # 3) gripper arm reaches out, grabs and opens the panel, holds it, releases
    {"joint": "shoulder_1", "unit": "deg",
     "keys": [[0, 0], [5, 0], [6.5, -62], [8, -34], [13, -34], [15, 0], [16, 0]]},
    {"joint": "elbow_1", "unit": "deg",
     "keys": [[0, 0], [5, 0], [6.5, 85], [8, 60], [13, 60], [15, 0], [16, 0]]},
    {"joint": "wrist_1", "unit": "deg",
     "keys": [[0, 0], [6.5, 0], [7.5, 28], [13, 28], [15, 0], [16, 0]]},
    {"joint": "gripper_1", "unit": "deg",
     "keys": [[0, 0], [6, 0], [6.6, 36], [8, 36], [12.5, 36], [13, 4], [15, 0], [16, 0]]},

    # 4) camera arm swings in to inspect the opened interior, then withdraws
    {"joint": "shoulder_2", "unit": "deg",
     "keys": [[0, 0], [8, 0], [9.5, -52], [12, -52], [14.5, 0], [16, 0]]},
    {"joint": "elbow_2", "unit": "deg",
     "keys": [[0, 0], [8, 0], [9.5, 72], [12, 72], [14.5, 0], [16, 0]]},
    {"joint": "wrist_2", "unit": "deg",
     "keys": [[0, 0], [9.5, 0], [10.5, -42], [12, -42], [14.5, 0], [16, 0]]},

    # 5) wrench arm drives down and the socket spins to tighten one fastener
    {"joint": "shoulder_3", "unit": "deg",
     "keys": [[0, 0], [9, 0], [10.5, -56], [13.5, -56], [15, 0], [16, 0]]},
    {"joint": "elbow_3", "unit": "deg",
     "keys": [[0, 0], [9, 0], [10.5, 92], [13.5, 92], [15, 0], [16, 0]]},
    {"joint": "wrist_3", "unit": "deg",
     "keys": [[0, 0], [10.5, 0], [11.5, 18], [13.5, 18], [15, 0], [16, 0]]},
    {"joint": "wrench_spin_3", "unit": "deg",
     "keys": [[0, 0], [11, 0], [13.3, 360], [16, 360]]},
]}